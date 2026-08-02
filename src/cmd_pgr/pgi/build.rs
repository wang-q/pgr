//! `pgr pgi build` — build a .pgi index from FASTA or 2bit.

use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};

/// Build the clap subcommand for build.
pub fn make_subcommand() -> Command {
    Command::new("build")
        .about("Builds a pgr genome index (.pgi) from FASTA or 2bit")
        .after_help(
            r###"
Builds a syncmer-sparse sorted k-mer index. Each syncmer position (default
(12,8) canonical, matching FastGA GIX) seeds a k-mer of length -k on both
strands (unless --no-rev). The result is a single binary file, usable for
`pgr dist pgi` (merge distance), `pgr pgi to-hv` (hypervector projection) and
future FastGA-style seed discovery.

Notes:
* Input FASTA may be plain or gzipped; .2bit input is supported (fastest).
* K-mers containing non-ACGT bases (e.g. N) are skipped.
* --mask skips soft-masked regions (lowercase FASTA bases / 2bit mask blocks),
  matching FastGA -M: masked repeats and low-complexity regions produce no seeds.

Examples:
1. Build from FASTA:
   pgr pgi build genome.fa.gz -o genome.pgi
2. Build from 2bit:
   pgr pgi build genome.2bit -o genome.pgi
"###,
        )
        .arg(
            Arg::new("infile")
                .index(1)
                .required(true)
                .help("FASTA or .2bit genome file"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg_required())
        .arg(
            Arg::new("k")
                .short('k')
                .long("kmer")
                .default_value("40")
                .value_parser(value_parser!(usize))
                .help("k-mer size (must be <= 64)"),
        )
        .arg(
            Arg::new("smer")
                .long("smer")
                .default_value("8")
                .value_parser(value_parser!(usize))
                .help("Syncmer s-mer length"),
        )
        .arg(
            Arg::new("window")
                .long("window")
                .default_value("5")
                .value_parser(value_parser!(usize))
                .help("Syncmer window (s-mers); span = smer + window - 1"),
        )
        .arg(
            Arg::new("no_rev")
                .long("no-rev")
                .action(ArgAction::SetTrue)
                .help("Index the forward strand only (default: both)"),
        )
        .arg(
            Arg::new("mask")
                .long("mask")
                .action(ArgAction::SetTrue)
                .help("Skip soft-masked regions (lowercase FASTA / 2bit mask blocks)"),
        )
}
/// Execute the build command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = args.get_one::<String>("outfile").unwrap();
    let k = *args.get_one::<usize>("k").unwrap();
    let smer = *args.get_one::<usize>("smer").unwrap();
    let window = *args.get_one::<usize>("window").unwrap();
    let no_rev = args.get_flag("no_rev");
    let mask = args.get_flag("mask");

    let idx = pgr::libs::pgi::build::build_from_path(infile, k, smer, window, no_rev, mask)?;
    let mut writer = pgr::writer(outfile)?;
    idx.write(&mut writer)?;
    log::info!(
        "wrote {} unique k-mers / {} positions (k={}, syncmer {smer}/{window}) to {}",
        idx.n_unique(),
        idx.n_positions(),
        k,
        outfile
    );
    Ok(())
}
