use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for sample.
pub fn make_subcommand() -> Command {
    Command::new("sample")
        .about("Subsamples reads to a target base count")
        .after_help(
            r###"
This command downsamples reads so the output contains approximately the
requested number of bases, preserving input order. Selection is deterministic
for a given seed and matches BBTools `reformat.sh` with `samplebasestarget`
and a fixed `sampleseed`.

Notes:
* Requires a file input: two passes are made over the data
* Supports both plain text and gzipped (.gz) files
* The default seed (1) reproduces BBTools sampleseed=1 selection

Examples:
1. Keep about 1 million bases:
   pgr fq sample reads.fq -o out.fq --bases 1000000

2. Reproduce a BBTools run with a fixed seed:
   pgr fq sample reads.fq.gz -o out.fq --bases 1000000 --seed 42
"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg_with_numargs(
            "Input FASTQ file",
            1..=1,
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("bases")
                .long("bases")
                .num_args(1)
                .value_parser(value_parser!(i64))
                .help("Target number of output bases"),
        )
        .arg(
            Arg::new("seed")
                .long("seed")
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(u64))
                .help("Random seed for deterministic selection"),
        )
}

/// Execute the sample command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infiles").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let bases = *args.get_one::<i64>("bases").unwrap();
    let seed = *args.get_one::<u64>("seed").unwrap();
    if bases < 0 {
        anyhow::bail!("--bases must be non-negative");
    }
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, [infile.as_str()])?;
    let mut out =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    pgr::libs::fq::sample::sample(infile, &mut out, bases, seed)?;
    out.flush()?;
    Ok(())
}
