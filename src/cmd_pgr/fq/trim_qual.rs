use anyhow::{bail, Context};
use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};
use pgr::libs::fq::trim::{Method, QualityBase, TrimOptions};
use std::io::Write;

/// Build the clap subcommand for trim-qual.
pub fn make_subcommand() -> Command {
    Command::new("trim-qual")
        .about("Trims reads by quality score")
        .after_help(
            r###"
This command trims low-quality bases from read ends using a sliding window
(default) or the Mott cumulative-quality algorithm.

Notes:
* Quality trimming only; adapters are not removed
* Supports single-end (1 file) and paired-end (2 files) FASTQ input
* For paired-end input, omit --outfile-2 to write interleaved output
* With --outfile-single, surviving reads whose mate failed are written there
* Reads shorter than the length threshold after trimming are discarded
* Quality encoding is auto-detected (33/64) by default; override with --quality-base
* Supports both plain text and gzipped (.gz) files

Examples:
1. Single-end trimming:
   pgr fq trim-q in.fq -o out.fq

2. Paired-end with separate outputs and singles:
   pgr fq trim-q R1.fq R2.fq -o r1.fq --outfile-2 r2.fq --outfile-single s.fq

3. Paired-end interleaved output:
   pgr fq trim-q R1.fq R2.fq -o interleaved.fq
"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg_with_numargs(
            "Input FASTQ file(s): 1 single-end or 2 paired-end",
            1..=2,
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("outfile_2")
                .long("outfile-2")
                .num_args(1)
                .help("R2 output file (paired-end; omit for interleaved output)"),
        )
        .arg(
            Arg::new("outfile_single")
                .long("outfile-single")
                .num_args(1)
                .help("Output file for surviving single-end reads"),
        )
        .arg(
            Arg::new("qual")
                .long("qual-threshold")
                .short('q')
                .num_args(1)
                .default_value("20")
                .value_parser(value_parser!(f64))
                .help("Quality threshold"),
        )
        .arg(
            Arg::new("length")
                .long("length-threshold")
                .short('l')
                .num_args(1)
                .default_value("20")
                .value_parser(value_parser!(usize))
                .help("Minimum kept length; shorter reads are discarded"),
        )
        .arg(
            Arg::new("method")
                .long("method")
                .num_args(1)
                .default_value("sliding")
                .value_parser(["sliding", "mott"])
                .help("Trimming algorithm"),
        )
        .arg(
            Arg::new("no_fiveprime")
                .long("no-fiveprime")
                .action(ArgAction::SetTrue)
                .help("Disable 5' trimming"),
        )
        .arg(
            Arg::new("quality_base")
                .long("quality-base")
                .num_args(1)
                .default_value("auto")
                .value_parser(["33", "64", "auto"])
                .help("Input quality encoding (auto-detected by default)"),
        )
        .arg(
            Arg::new("polyg_right")
                .long("polyg-right")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Trim 3' poly-G tails of at least this length (0 disables)"),
        )
}

/// Execute the trim-q command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let outfile_2 = args.get_one::<String>("outfile_2").map(String::as_str);
    let outfile_single = args.get_one::<String>("outfile_single").map(String::as_str);

    if infiles.len() == 1 && outfile_2.is_some() {
        bail!("--outfile-2 requires two input files (paired-end)");
    }
    for out in [outfile_2, outfile_single].into_iter().flatten() {
        if out == "stdout" {
            bail!("--outfile-2/--outfile-single must be file paths, not stdout");
        }
    }
    for out in [Some(outfile), outfile_2, outfile_single]
        .into_iter()
        .flatten()
    {
        crate::cmd_pgr::args::ensure_outfile_distinct(out, infiles.iter().map(String::as_str))?;
    }
    let outputs = [Some(outfile), outfile_2, outfile_single]
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
    for (i, a) in outputs.iter().enumerate() {
        for b in &outputs[i + 1..] {
            if *a != "stdout" && *b != "stdout" && pgr::libs::io::same_path(a, b) {
                bail!("output files must be distinct: {} and {}", a, b);
            }
        }
    }

    let method = if args.get_one::<String>("method").unwrap() == "mott" {
        Method::Mott
    } else {
        Method::Sliding
    };
    let quality_base = match args.get_one::<String>("quality_base").unwrap().as_str() {
        "33" => QualityBase::Phred33,
        "64" => QualityBase::Phred64,
        _ => QualityBase::Auto,
    };
    let opts = TrimOptions {
        qual_threshold: *args.get_one::<f64>("qual").unwrap(),
        length_threshold: *args.get_one::<usize>("length").unwrap(),
        method,
        no_fiveprime: args.get_flag("no_fiveprime"),
        quality_base,
        polyg_right: *args.get_one::<usize>("polyg_right").unwrap(),
    };

    if infiles.len() == 1 {
        let mut out = pgr::writer(outfile)
            .with_context(|| format!("Failed to open writer for {}", outfile))?;
        pgr::libs::fq::trim::run_single(&infiles[0], &mut out, &opts)?;
        out.flush()?;
    } else {
        let mut out1 = pgr::writer(outfile)
            .with_context(|| format!("Failed to open writer for {}", outfile))?;
        let mut out2 = outfile_2
            .map(|p| pgr::writer(p).with_context(|| format!("Failed to open writer for {}", p)))
            .transpose()?;
        let mut single = outfile_single
            .map(|p| pgr::writer(p).with_context(|| format!("Failed to open writer for {}", p)))
            .transpose()?;
        pgr::libs::fq::trim::run_paired(
            &infiles[0],
            &infiles[1],
            &mut out1,
            out2.as_mut(),
            single.as_mut(),
            &opts,
        )?;
        out1.flush()?;
        if let Some(w) = out2.as_mut() {
            w.flush()?;
        }
        if let Some(w) = single.as_mut() {
            w.flush()?;
        }
    }
    Ok(())
}
