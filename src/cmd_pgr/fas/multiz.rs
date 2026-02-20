use clap::*;
use std::io::Write;

pub fn make_subcommand() -> Command {
    Command::new("multiz")
        .about("Merge block FA files using multiz-like DP on reference")
        .after_help(
            r###"
Merge multiple block FA files in the shared reference coordinate system using a multiz-like banded DP.

Notes:
* Takes two or more .fas inputs that share a reference name.
* Automatically derives windows from reference coverage with radius padding.
* Supports core (intersection) and union modes on windows and species.

Examples:
1. Core mode merge with default radius:
   pgr fas multiz -r ref tests/fas/part1.fas tests/fas/part2.fas

2. Union mode with larger radius and minimum width:
   pgr fas multiz -r ref --mode union --radius 30 --min-width 1000 part1.fas part2.fas part3.fas

3. Write merged blocks to a file:
   pgr fas multiz -r ref part1.fas part2.fas -o merged.fas
"###,
        )
        .arg(
            Arg::new("ref")
                .short('r')
                .long("ref")
                .num_args(1)
                .required(true)
                .help("Reference sequence name present in all inputs"),
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(2..)
                .index(1)
                .help("Input block FA file(s) to merge"),
        )
        .arg(
            Arg::new("radius")
                .long("radius")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("30")
                .help("Banded DP radius around the reference diagonal"),
        )
        .arg(
            Arg::new("min_width")
                .long("min-width")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("1")
                .help("Minimum window width to consider for merging"),
        )
        .arg(
            Arg::new("mode")
                .long("mode")
                .num_args(1)
                .default_value("core")
                .help("Merge mode: core (strict intersection) or union"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let ref_name = args.get_one::<String>("ref").unwrap().to_string();
    let radius = *args.get_one::<usize>("radius").unwrap();
    let min_width = *args.get_one::<usize>("min_width").unwrap();
    let mode_str = args.get_one::<String>("mode").unwrap();

    let mode = match mode_str.as_str() {
        "core" => pgr::libs::fas_multiz::FasMultizMode::Core,
        "union" => pgr::libs::fas_multiz::FasMultizMode::Union,
        other => {
            return Err(anyhow::anyhow!(format!(
                "Invalid mode '{}', expected 'core' or 'union'",
                other
            )))
        }
    };

    let cfg = pgr::libs::fas_multiz::FasMultizConfig {
        ref_name: ref_name.clone(),
        radius,
        min_width,
        mode,
        match_score: 2,
        mismatch_score: -1,
        gap_score: -2,
    };

    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();

    let blocks = pgr::libs::fas_multiz::merge_fas_files_auto_windows(&ref_name, &infiles, &cfg)?;

    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap());

    for block in blocks {
        for entry in &block.entries {
            let range = entry.range();
            let seq = String::from_utf8(entry.seq().to_vec()).unwrap();
            writeln!(writer, ">{}", range)?;
            writeln!(writer, "{}", seq)?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

