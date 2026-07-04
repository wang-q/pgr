use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};

/// Build the clap subcommand for to-range.
pub fn make_subcommand() -> Command {
    Command::new("to-range")
        .about("Extracts coordinates from PSL as ranges (.rg)")
        .after_help(
            r###"
Extract alignment coordinates from PSL files and output in .rg format (chr:start-end).
This is useful for depth calculation with `spanr coverage`.

Notes:
* Coordinates are converted to 1-based inclusive.
* Supports strand-aware coordinate conversion (outputs positive strand coordinates).
* Outputs one range per alignment block.

Examples:
1. Extract query ranges:
   pgr psl to-range input.psl > query.rg

2. Extract target ranges:
   pgr psl to-range input.psl --target-coords > target.rg
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input PSL file. [stdin] for standard input",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("target_coords")
                .long("target-coords")
                .action(ArgAction::SetTrue)
                .help("Extract target coordinates instead of query"),
        )
        .arg(
            Arg::new("strict")
                .long("strict")
                .action(ArgAction::SetTrue)
                .help("Fail on parse errors instead of skipping malformed lines"),
        )
}

/// Execute the to-range command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let output = crate::cmd_pgr::args::get_outfile(args);
    let extract_target = args.get_flag("target_coords");
    let strict = args.get_flag("strict");

    let reader =
        pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;
    let mut writer =
        pgr::writer(output).with_context(|| format!("Failed to open writer for {}", output))?;

    pgr::libs::fmt::psl::to_ranges(reader, &mut writer, extract_target, strict)?;

    Ok(())
}
