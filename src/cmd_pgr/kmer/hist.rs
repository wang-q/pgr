use clap::{Arg, ArgMatches, Command};

/// Build the clap subcommand for hist.
pub fn make_subcommand() -> Command {
    Command::new("hist")
        .about("Builds a k-mer frequency histogram (.hist)")
        .after_help(
            r###"
Builds a k-mer frequency histogram from sequences or an existing table and
writes it in the FastK `.hist` binary layout (readable by Histex, KatGC,
and GenomeScope tooling).

Give either a sequence file (histogram is computed on the fly) or --table
to reuse an existing `.pkt` table. Bins are fixed 1..=32767; counts above
the top bin are folded into it, matching FastK semantics.

* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

Examples:
1. Histogram from reads:
   pgr kmer hist reads.fq.gz -k 21 -o reads.hist
2. Histogram from an existing table:
   pgr kmer hist -t reads.pkt -o reads.hist
"###,
        )
        .arg(
            Arg::new("infile")
                .num_args(1)
                .index(1)
                .required_unless_present("table")
                .help("Input FASTA/FASTQ file to process (unless --table is given)"),
        )
        .arg(super::profile::table_arg())
        .arg(super::profile::kmer_arg())
        .arg(crate::cmd_pgr::args::outfile_arg_required())
}

/// Execute the hist command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = args.get_one::<String>("outfile").unwrap();
    let table_path = args.get_one::<String>("table").map(String::as_str);
    let k = super::resolve_k(args.get_one::<usize>("kmer"), table_path)?;

    let table = if let Some(t) = table_path {
        pgr::libs::kmer::count::load(std::path::Path::new(t), k)?
    } else {
        let infile = args.get_one::<String>("infile").unwrap();
        let seqs = super::read_seqs(infile)?;
        pgr::libs::kmer::count::build_table(&seqs, k)?
    };
    let hist = pgr::libs::kmer::hist::from_table(&table);
    pgr::libs::kmer::hist::write(std::path::Path::new(outfile), &hist)?;
    log::info!(
        "==> Wrote histogram of {} distinct {}-mers to {}",
        table.keys.len(),
        k,
        outfile
    );
    Ok(())
}
