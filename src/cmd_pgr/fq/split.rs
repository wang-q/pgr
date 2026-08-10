use anyhow::{bail, Context};
use clap::{Arg, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for split.
pub fn make_subcommand() -> Command {
    Command::new("split")
        .about("Splits interleaved FASTQ into R1/R2/singles files")
        .after_help(
            r###"
This command splits an interleaved FASTQ file into paired-end R1/R2 outputs
and a singles file for unpaired reads. It is the inverse of `pgr fq interleave`
and matches BBTools `repair.sh` in `rp` mode. By default reads are paired by
position (every two records); `--repair` instead matches mates by read name
prefix, recovering disordered pairs and routing orphaned reads to singles.

Notes:
* Reads are processed in input order; headers and quality are preserved
* A trailing read without its mate is written to --outfile-single
* --repair buffers unpaired reads in memory (like repair.sh); the default
  positional mode is streaming
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

Examples:
1. Split into R1/R2 and singles:
   pgr fq split interleaved.fq -o r1.fq --outfile-2 r2.fq --outfile-single s.fq

2. Split into R1/R2 only:
   pgr fq split interleaved.fq.gz -o r1.fq --outfile-2 r2.fq
"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg_with_numargs(
            "Input interleaved FASTQ file",
            1..=1,
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("outfile_2")
                .long("outfile-2")
                .num_args(1)
                .required(true)
                .help("R2 output file (required)"),
        )
        .arg(
            Arg::new("outfile_single")
                .long("outfile-single")
                .num_args(1)
                .help("Output file for unpaired reads"),
        )
        .arg(
            Arg::new("repair")
                .long("repair")
                .action(clap::ArgAction::SetTrue)
                .help("Pair mates by read-name prefix instead of position"),
        )
}

/// Execute the split command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infiles").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let outfile_2 = args.get_one::<String>("outfile_2").unwrap();
    let outfile_single = args.get_one::<String>("outfile_single");

    if outfile_2 == "stdout" {
        bail!("--outfile-2 must be a file path, not stdout");
    }
    if outfile_single == Some(&"stdout".to_string()) {
        bail!("--outfile-single must be a file path, not stdout");
    }
    for out in [
        Some(outfile),
        Some(outfile_2.as_str()),
        outfile_single.map(String::as_str),
    ]
    .into_iter()
    .flatten()
    {
        crate::cmd_pgr::args::ensure_outfile_distinct(out, [infile.as_str()])?;
    }

    let mut out1 =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    let mut out2 = pgr::writer(outfile_2)
        .with_context(|| format!("Failed to open writer for {}", outfile_2))?;
    let mut single = outfile_single
        .map(|p| pgr::writer(p).with_context(|| format!("Failed to open writer for {}", p)))
        .transpose()?;
    if args.get_flag("repair") {
        pgr::libs::fq::split::split_repair(infile, &mut out1, &mut out2, single.as_mut())?;
    } else {
        pgr::libs::fq::split::split(infile, &mut out1, &mut out2, single.as_mut())?;
    }
    out1.flush()?;
    out2.flush()?;
    if let Some(w) = single.as_mut() {
        w.flush()?;
    }
    Ok(())
}
