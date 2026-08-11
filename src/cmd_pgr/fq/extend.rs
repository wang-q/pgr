use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use pgr::libs::asm::tadpole::{run, TadpoleOptions};
use std::io::Write;

/// Build the clap subcommand for extend.
pub fn make_subcommand() -> Command {
    Command::new("extend")
        .about("Extends reads along the k-mer graph (tadpole-compatible)")
        .after_help(
            r###"
This command extends reads in both directions along the k-mer graph, stopping
at junctions and dead ends, reproducing the BBTools `tadpole.sh mode=extend`
behavior (k>31 uses the Tadpole2 long-k-mer path). Unlike `fq ecc`, extend
mode does not run k-mer error correction.

Notes:
* Input is 1 interleaved FASTQ file or 2 paired files (R1, R2)
* Extended bases get BBTools' fake quality (phred 30); other qualities are
  canonicalized like BBTools (phred round-trip)
* `--extend-rollback` trims random trailing bases of partial extensions
  (deterministic, based on the read's input index)
* Processing is ordered and deterministic (equivalent to `threads=1`)
* Supports both plain text and gzipped (.gz) files

Examples:
1. Extend by 20 bp each side with k=62 (anchr read-extension step):
   pgr fq extend in.fq.gz -o extended.fq.gz --kmer 62 --el 20 --er 20

2. Extend only to the right:
   pgr fq extend in.fq.gz -o out.fq --el 0 --er 50
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
                .value_parser(clap::builder::RangedU64ValueParser::<usize>::new().range(1..))
                .help("K-mer length"),
        )
        .arg(
            Arg::new("el")
                .long("el")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Extend to the left by at most this many bases"),
        )
        .arg(
            Arg::new("er")
                .long("er")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Extend to the right by at most this many bases"),
        )
        .arg(
            Arg::new("min_prob")
                .long("min-prob")
                .num_args(1)
                .value_parser(value_parser!(f32))
                .help("Ignore k-mers below this error-free probability"),
        )
        .arg(
            Arg::new("extend_rollback")
                .long("extend-rollback")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Trim up to this many bases of partial extensions"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .num_args(1)
                .default_value("auto")
                .help("Accepted for tadpole.sh compatibility; ignored (deterministic single-pass)"),
        )
}

/// Execute the extend command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut opts = TadpoleOptions {
        k: *args.get_one::<usize>("kmer").unwrap(),
        extend_left: args.get_one::<usize>("el").copied().unwrap_or(100),
        extend_right: args.get_one::<usize>("er").copied().unwrap_or(100),
        ..Default::default()
    };
    // Validate the thread-count value; processing stays deterministic
    // single-pass (see the design notes), so the result is not used.
    crate::cmd_pgr::args::parse_parallel_auto(args.get_one::<String>("parallel").unwrap())?;
    if let Some(x) = args.get_one::<f32>("min_prob") {
        opts.min_prob = *x;
    }
    if let Some(x) = args.get_one::<usize>("extend_rollback") {
        opts.extension_rollback = *x;
    }

    let mut out = pgr::libs::io::writer(outfile)
        .with_context(|| format!("failed to open output {outfile}"))?;
    let stats = run(&infiles, &mut out, &opts)?;
    out.flush()?;
    eprintln!(
        "Reads in: {}  Reads extended: {}  Bases extended: {}  Reads discarded: {}",
        stats.reads_in, stats.reads_extended, stats.bases_extended, stats.reads_discarded
    );
    Ok(())
}
