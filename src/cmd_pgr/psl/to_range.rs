use clap::*;
use pgr::libs::fmt::psl::Psl;
use std::io::{BufRead, Write};
use std::str::FromStr;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("to-range")
        .about("Extract coordinates from PSL as ranges (.rg)")
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
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Input PSL file. [stdin] for standard input"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("target_coords")
                .long("target-coords")
                .short('t')
                .action(ArgAction::SetTrue)
                .help("Extract target coordinates instead of query"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;
    let infile = args.get_one::<String>("infile").unwrap();
    let reader = pgr::reader(infile)?;
    let extract_target = args.get_flag("target_coords");

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        // Skip header lines
        if line.starts_with("psLayout") || line.starts_with("match") || line.starts_with("------") {
            continue;
        }

        let psl = match Psl::from_str(&line) {
            Ok(p) => p,
            Err(_) => continue,
        };

        for range in pgr::libs::fmt::psl::psl_block_ranges(&psl, extract_target) {
            writer.write_all(range.as_bytes())?;
            writer.write_all(b"\n")?;
        }
    }

    Ok(())
}
