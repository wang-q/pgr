use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use pgr::libs::fq::clump::{clump, ClumpOptions};
use std::io::Write;

/// Build the clap subcommand for clump.
pub fn make_subcommand() -> Command {
    Command::new("clump")
        .about("Sorts reads by k-mer signature (clumpify-compatible)")
        .after_help(
            r###"
This command sorts interleaved paired reads by the pivot k-mer of R1,
reproducing the BBTools `clumpify.sh` default output order byte for byte.
The sorting clusters reads that share k-mers, which speeds up the k-mer
steps that follow in a read-cleaning pipeline.

Notes:
* Paired input must be interleaved; mates are kept together
* Deterministic for a given k-mer size and seed
* --dedupe removes whole-pair duplicates (R1 and R2 both matching within
  --dupesubs substitutions), keeping the higher-quality copy
* --mem sets the in-memory sort budget (KMG, default 2g); data estimated to
  exceed it is sorted via external hash buckets (--buckets to override the
  bucket count)
* --sort-mode forces the path: auto (default, by memory budget), global
  (always in-memory), or bucket (always external); specifying --buckets
  implies bucket mode
* --parallel caps the parallel worker pool (default: logical CPU count)
* Supports both plain text and gzipped (.gz) files

Examples:
1. Sort reads with the BBTools-compatible defaults:
   pgr fq clump in.fq.gz -o clumped.fq

2. Reproduce a BBTools run with a different seed:
   pgr fq clump in.fq.gz -o out.fq --seed 2

3. Remove exact duplicate pairs:
   pgr fq clump R1.fq.gz R2.fq.gz -o out.fq --dedupe --dupesubs 0

4. Bound memory to 1 GiB (external bucket path for larger data):
   pgr fq clump R1.fq.gz R2.fq.gz -o out.fq --mem 1g

5. Force the external bucket path:
   pgr fq clump R1.fq.gz R2.fq.gz -o out.fq --sort-mode bucket
"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg_with_numargs(
            "Input FASTQ file(s): 1 interleaved or 2 paired (R1, R2)",
            1..=2,
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("kmer")
                .long("kmer")
                .short('k')
                .num_args(1)
                .default_value("31")
                .value_parser(value_parser!(usize))
                .help("K-mer size (2..=31)"),
        )
        .arg(
            Arg::new("seed")
                .long("seed")
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(u64))
                .help("Comparator seed"),
        )
        .arg(
            Arg::new("dedupe")
                .long("dedupe")
                .action(clap::ArgAction::SetTrue)
                .help("Remove duplicate read pairs"),
        )
        .arg(
            Arg::new("dupesubs")
                .long("dupesubs")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Maximum substitutions allowed in a duplicate"),
        )
        .arg(
            Arg::new("mem")
                .long("mem")
                .num_args(1)
                .default_value("2g")
                .help("In-memory sort budget (KMG; default 2g)"),
        )
        .arg(
            Arg::new("buckets")
                .long("buckets")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("External-path hash bucket count"),
        )
        .arg(
            Arg::new("sort_mode")
                .long("sort-mode")
                .num_args(1)
                .default_value("auto")
                .value_parser(["auto", "global", "bucket"])
                .help("Sorting path: auto (default), global, or bucket"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .num_args(1)
                .default_value("auto")
                .help("Worker threads (default: logical CPU count)"),
        )
}

/// Execute the clump command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let k = *args.get_one::<usize>("kmer").unwrap();
    let seed = *args.get_one::<u64>("seed").unwrap();
    let dedupe = args.get_flag("dedupe");
    let dupesubs = *args.get_one::<usize>("dupesubs").unwrap();
    let mem = Some(pgr::libs::sys::parse_mem_size(
        args.get_one::<String>("mem").unwrap(),
    )?);
    let buckets = args.get_one::<usize>("buckets").copied();
    let mode = match args.get_one::<String>("sort_mode").unwrap().as_str() {
        "global" => pgr::libs::fq::clump::SortMode::Global,
        "bucket" => pgr::libs::fq::clump::SortMode::Bucket,
        _ => pgr::libs::fq::clump::SortMode::Auto,
    };
    // Specifying a bucket count implies the external path.
    let mode = if mode == pgr::libs::fq::clump::SortMode::Auto && buckets.is_some() {
        pgr::libs::fq::clump::SortMode::Bucket
    } else {
        mode
    };
    let parallel =
        crate::cmd_pgr::args::parse_parallel_auto(args.get_one::<String>("parallel").unwrap())?;
    if !(2..=31).contains(&k) {
        anyhow::bail!("--kmer must be in 2..=31, got {}", k);
    }
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, infiles.iter().map(String::as_str))?;
    let mut out =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    let opts = ClumpOptions {
        k,
        seed,
        dedupe,
        dupesubs,
        mem,
        buckets,
        mode,
        parallel,
    };
    clump(&infiles, &mut out, &opts)?;
    out.flush()?;
    Ok(())
}
