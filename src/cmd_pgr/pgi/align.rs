//! `pgr pgi align` — two-index merge alignment into PSL blocks.

use clap::{value_parser, Arg, ArgMatches, Command};

/// Build the clap subcommand for align.
pub fn make_subcommand() -> Command {
    Command::new("align")
        .about("Aligns two .pgi indexes into PSL blocks")
        .after_help(
            r###"
Merges two sorted .pgi k-mer streams, chains the shared seeds in
anti-diagonal space, and emits one PSL block per chain. Block-level output is
meant to be chained by `pgr psl to_chain` / `pgr pl chainnet`.
With --ref-seq/--query-seq, each chain is refined by a banded local alignment
into a scored PSL record with real blocks (chains longer than 30 kb stay as
plain blocks).

Notes:
* Both indexes must use identical sampling parameters (k, syncmer, window).
* K-mers occurring more than --freq times on either side are skipped.
* Chains shorter than --min-span on either axis are dropped.
* --ref-seq/--query-seq accept FASTA (.fa/.fa.gz) or .2bit files whose
  sequences correspond by index to the .pgi contig order.

Examples:
1. Align two indexes:
   pgr pgi align ref.pgi query.pgi -o out.psl
2. Tune seed filtering and chaining:
   pgr pgi align ref.pgi query.pgi -f 20 -c 100 -s 2000 --band 64 -o out.psl
3. Refine chains with the source sequences:
   pgr pgi align ref.pgi query.pgi --ref-seq ref.fa --query-seq query.fa -o out.psl
4. Stitch chains across small insertions:
   pgr pgi align ref.pgi query.pgi --merge-gap 10000 -o out.psl
5. Partial seeds are experimental and degrade results; exact k-mers default:
   pgr pgi align ref.pgi query.pgi --min-shared 30 -o out.psl
"###,
        )
        .arg(
            Arg::new("ref_idx")
                .index(1)
                .required(true)
                .help("Reference .pgi index"),
        )
        .arg(
            Arg::new("query_idx")
                .index(2)
                .required(true)
                .help("Query .pgi index"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg_required())
        .arg(
            Arg::new("freq")
                .short('f')
                .long("freq")
                .default_value("10")
                .value_parser(value_parser!(u32))
                .help("Maximum k-mer frequency on either side to keep as seed"),
        )
        .arg(
            Arg::new("min_span")
                .short('c')
                .long("min-span")
                .default_value("85")
                .value_parser(value_parser!(u32))
                .help("Minimum per-axis seed span (bp) for a chain"),
        )
        .arg(
            Arg::new("max_gap")
                .short('s')
                .long("max-gap")
                .default_value("1000")
                .value_parser(value_parser!(u32))
                .help("Maximum bp gap between consecutive seeds in a chain"),
        )
        .arg(
            Arg::new("band")
                .long("band")
                .default_value("128")
                .value_parser(value_parser!(u32))
                .help("Diagonal band half-width (bp) around the chain mean"),
        )
        .arg(
            Arg::new("merge_gap")
                .long("merge-gap")
                .default_value("5000")
                .value_parser(value_parser!(u32))
                .help("Maximum gap (bp) between adjacent colinear chains to merge"),
        )
        .arg(
            Arg::new("min_shared")
                .long("min-shared")
                .value_parser(value_parser!(usize))
                .help("Minimum shared seed length (bp); default = k for greedy, k/2 for tube"),
        )
        .arg(
            Arg::new("workflow")
                .long("workflow")
                .default_value("greedy")
                .value_parser(["greedy", "tube"])
                .help("Chaining workflow: greedy chains (default) or FastGA tubes"),
        )
        .arg(
            Arg::new("ref_seq")
                .long("ref-seq")
                .help("Reference sequence file (FASTA or .2bit) for chain refinement"),
        )
        .arg(
            Arg::new("query_seq")
                .long("query-seq")
                .help("Query sequence file (FASTA or .2bit) for chain refinement"),
        )
}

/// Execute the align command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let ref_idx = args.get_one::<String>("ref_idx").unwrap();
    let query_idx = args.get_one::<String>("query_idx").unwrap();
    let outfile = args.get_one::<String>("outfile").unwrap();
    let params = pgr::libs::pgi::align::AlignParams {
        freq: *args.get_one::<u32>("freq").unwrap(),
        min_span: *args.get_one::<u32>("min_span").unwrap(),
        max_gap: *args.get_one::<u32>("max_gap").unwrap(),
        band: *args.get_one::<u32>("band").unwrap(),
        merge_gap: *args.get_one::<u32>("merge_gap").unwrap(),
        min_shared: args.get_one::<usize>("min_shared").copied(),
        workflow: match args.get_one::<String>("workflow").unwrap().as_str() {
            "tube" => pgr::libs::pgi::align::Workflow::Tube,
            _ => pgr::libs::pgi::align::Workflow::Greedy,
        },
    };

    let mut r1 = pgr::reader(ref_idx)?;
    let a = pgr::libs::pgi::PgiIndex::read(&mut r1)?;
    let mut r2 = pgr::reader(query_idx)?;
    let b = pgr::libs::pgi::PgiIndex::read(&mut r2)?;

    let psls = match (
        args.get_one::<String>("ref_seq"),
        args.get_one::<String>("query_seq"),
    ) {
        (Some(rp), Some(qp)) => {
            let rs = read_seqs(rp)?;
            let qs = read_seqs(qp)?;
            pgr::libs::pgi::align::align_to_psl_ext(&a, &b, &params, &rs, &qs)?
        }
        _ => pgr::libs::pgi::align::align_to_psl(&a, &b, &params)?,
    };
    let mut writer = pgr::writer(outfile)?;
    for p in &psls {
        p.write_to(&mut writer)?;
    }
    log::info!(
        "wrote {} PSL blocks (freq={}, min-span={}, max-gap={}, band={}) to {}",
        psls.len(),
        params.freq,
        params.min_span,
        params.max_gap,
        params.band,
        outfile
    );
    Ok(())
}

/// Read all sequences from a FASTA (plain or gzipped) or .2bit file.
fn read_seqs(path: &str) -> anyhow::Result<Vec<(String, Vec<u8>)>> {
    let is_2bit = std::path::Path::new(path)
        .extension()
        .and_then(|e| e.to_str())
        == Some("2bit");
    if is_2bit {
        pgr::libs::pgi::build::read_2bit(path)
    } else {
        pgr::libs::pgi::build::read_fasta(path)
    }
}
