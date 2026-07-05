use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use std::io::Write;
/// Build the clap subcommand for multiz.
pub fn make_subcommand() -> Command {
    Command::new("multiz")
        .about("Merges block FA files using multiz-like DP on reference")
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
            Arg::new("ref_name")
                .short('r')
                .long("ref-name")
                .num_args(1)
                .required(true)
                .help("Reference sequence name present in all inputs"),
        )
        .arg(crate::cmd_pgr::args::infiles_arg_with_numargs(
            "Input block FA file(s) to merge",
            2..,
        ))
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
        .arg(crate::cmd_pgr::args::mode_arg(
            "core",
            &["core", "union"],
            "Merge mode: core (strict intersection) or union",
        ))
        .arg(crate::cmd_pgr::args::score_scheme_arg())
        .arg(crate::cmd_pgr::args::gap_model_arg(
            "medium",
            &["constant", "medium", "loose"],
            "Gap model: constant, medium, or loose",
        ))
        .arg(crate::cmd_pgr::args::align_gap_open_arg())
        .arg(crate::cmd_pgr::args::align_gap_extend_arg())
        .arg(
            Arg::new("match_score")
                .long("match-score")
                .num_args(1)
                .default_value("2")
                .value_parser(value_parser!(i32))
                .help("Match score for scoring matrix"),
        )
        .arg(
            Arg::new("mismatch_score")
                .long("mismatch-score")
                .num_args(1)
                .default_value("-1")
                .value_parser(value_parser!(i32))
                .help("Mismatch penalty for scoring matrix"),
        )
        .arg(
            Arg::new("gap_score")
                .long("gap-score")
                .num_args(1)
                .default_value("-2")
                .value_parser(value_parser!(i32))
                .help("Gap penalty for scoring matrix"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}
/// Execute the multiz command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let ref_name = args.get_one::<String>("ref_name").unwrap().to_string();
    let radius = *args.get_one::<usize>("radius").unwrap();
    let min_width = *args.get_one::<usize>("min_width").unwrap();
    let mode_str = args.get_one::<String>("mode").unwrap();
    let gap_model_str = args.get_one::<String>("gap_model").unwrap();
    let score_matrix = args.get_one::<String>("score_scheme").cloned();
    let gap_open = args.get_one::<i32>("align_gap_open").copied();
    let gap_extend = args.get_one::<i32>("align_gap_extend").copied();

    let mode = match mode_str.as_str() {
        "core" => pgr::libs::fas_multiz::FasMultizMode::Core,
        "union" => pgr::libs::fas_multiz::FasMultizMode::Union,
        _ => anyhow::bail!("unknown mode: {}", mode_str),
    };

    let gap_model = match gap_model_str.as_str() {
        "constant" => pgr::libs::fas_multiz::FasMultizGapModel::Constant,
        "medium" => pgr::libs::fas_multiz::FasMultizGapModel::Medium,
        "loose" => pgr::libs::fas_multiz::FasMultizGapModel::Loose,
        _ => anyhow::bail!("unknown gap_model: {}", gap_model_str),
    };

    let match_score = *args.get_one::<i32>("match_score").unwrap();
    let mismatch_score = *args.get_one::<i32>("mismatch_score").unwrap();
    let gap_score = *args.get_one::<i32>("gap_score").unwrap();

    let cfg = pgr::libs::fas_multiz::FasMultizConfig {
        ref_name: ref_name.clone(),
        radius,
        min_width,
        mode,
        match_score,
        mismatch_score,
        gap_score,
        gap_model,
        gap_open,
        gap_extend,
        score_matrix,
    };

    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();

    let blocks = pgr::libs::fas_multiz::merge_fas_files_auto_windows(&ref_name, &infiles, &cfg)?;

    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    for block in blocks {
        for entry in &block.entries {
            let range = entry.range();
            let seq = String::from_utf8(entry.seq().to_vec())?;
            writeln!(writer, ">{}", range)?;
            writeln!(writer, "{}", seq)?;
        }
        writeln!(writer)?;
    }

    writer.flush()?;
    Ok(())
}
