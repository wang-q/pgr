use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use pgr::libs::fq::tadpole::{run, TadpoleOptions};
use std::io::Write;

/// Build the clap subcommand for ec-kmer.
pub fn make_subcommand() -> Command {
    Command::new("ec-kmer")
        .about("Error-corrects reads by k-mer reassembly (tadpole-compatible)")
        .after_help(
            r###"
This command error-corrects reads through the k-mer graph (reassemble mode),
reproducing the BBTools `tadpole.sh ecc` behavior: k-mers are counted with a
quality gate (`--min-prob`), per-read errors are detected from k-mer depth
transitions and corrected by local reassembly, and reads can be discarded
with the `tossjunk` / `tossdepth` / `tossuncorrectable` flags.

Notes:
* Input is 1 interleaved FASTQ file or 2 paired files (R1, R2)
* Qualities are canonicalized like BBTools (phred round-trip), so output
  quality may differ from the input at N/low-quality positions
* Processing is ordered and deterministic (equivalent to `threads=1`)
* Supports both plain text and gzipped (.gz) files

Examples:
1. Error-correct with tadpole defaults (anchr merge phase 3):
   pgr fq ec-kmer in.fq.gz -o ecct.fq.gz --toss-junk --toss-depth 2 \
       --toss-uncorrectable

2. Only correct, keep everything:
   pgr fq ec-kmer R1.fq R2.fq -o corrected.fq --kmer 31
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
            Arg::new("min_prob")
                .long("min-prob")
                .num_args(1)
                .value_parser(value_parser!(f32))
                .help("Ignore k-mers below this error-free probability"),
        )
        .arg(
            Arg::new("toss_junk")
                .long("toss-junk")
                .action(clap::ArgAction::SetTrue)
                .help("Discard reads that cannot be used for assembly"),
        )
        .arg(
            Arg::new("toss_depth")
                .long("toss-depth")
                .num_args(1)
                .value_parser(value_parser!(i64))
                .help("Discard reads with k-mers at or below this depth"),
        )
        .arg(
            Arg::new("toss_uncorrectable")
                .long("toss-uncorrectable")
                .action(clap::ArgAction::SetTrue)
                .help("Discard reads with uncorrectable errors"),
        )
        .arg(
            Arg::new("low_depth_fraction")
                .long("low-depth-fraction")
                .num_args(1)
                .value_parser(value_parser!(f32))
                .help("Minimum low-depth k-mer fraction to discard a read"),
        )
        .arg(
            Arg::new("require_both_bad")
                .long("require-both-bad")
                .action(clap::ArgAction::SetTrue)
                .help("Only discard a pair if both reads fail"),
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

/// Execute the ec-kmer command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut opts = TadpoleOptions {
        k: *args.get_one::<usize>("kmer").unwrap(),
        ecc: true,
        ..Default::default()
    };
    // Validate the thread-count value; processing stays deterministic
    // single-pass (see the design notes), so the result is not used.
    crate::cmd_pgr::args::parse_parallel_auto(args.get_one::<String>("parallel").unwrap())?;
    if let Some(x) = args.get_one::<f32>("min_prob") {
        opts.min_prob = *x;
    }
    opts.toss_junk = args.get_flag("toss_junk");
    if let Some(x) = args.get_one::<i64>("toss_depth") {
        opts.toss_depth = *x;
    }
    opts.toss_uncorrectable = args.get_flag("toss_uncorrectable");
    if let Some(x) = args.get_one::<f32>("low_depth_fraction") {
        opts.low_depth_discard_fraction = *x;
    }
    opts.require_both_bad = args.get_flag("require_both_bad");

    let mut out = pgr::libs::io::writer(outfile)
        .with_context(|| format!("failed to open output {outfile}"))?;
    let stats = run(&infiles, &mut out, &opts)?;
    out.flush()?;
    eprintln!(
        "Reads in: {}  Reads detected: {}  Bases detected: {}  Bases corrected: {}  Reads corrected: {}  Reads discarded: {}  Rollbacks: {}",
        stats.reads_in, stats.reads_detected, stats.bases_detected, stats.bases_corrected, stats.reads_corrected, stats.reads_discarded, stats.rollbacks
    );
    Ok(())
}
