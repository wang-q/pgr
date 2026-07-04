use clap::{ArgMatches, Command};
use std::io::{BufRead, Write};
use std::str::FromStr;

use pgr::libs::fmt::psl::Psl;

/// Build the clap subcommand for lift.
pub fn make_subcommand() -> Command {
    Command::new("lift")
        .about("Lifts PSL coordinates from fragment alignments")
        .after_help(
            r###"
Lifts PSL coordinates from query/target fragments to genomic coordinates.

Notes:
* The query or target name must be in the format `chr:start-end`.
* The coordinates in the name are 1-based, inclusive.
* Requires a chromosome sizes file for correct negative strand lifting.

Examples:
1. Lift query coordinates:
   pgr psl lift input.psl --q-sizes chrom.sizes > output.psl

2. Lift target coordinates:
   pgr psl lift input.psl --t-sizes chrom.sizes > output.psl

3. Lift both:
   pgr psl lift input.psl --q-sizes q.sizes --t-sizes t.sizes > output.psl
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input PSL file. [stdin] for standard input",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(crate::cmd_pgr::args::q_sizes_arg())
        .arg(crate::cmd_pgr::args::t_sizes_arg())
}

/// Execute the lift command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;
    let infile = args.get_one::<String>("infile").unwrap();
    let reader = pgr::reader(infile)?;

    let q_sizes_file = args.get_one::<String>("q_sizes").map(|s| s.as_str());
    let t_sizes_file = args.get_one::<String>("t_sizes").map(|s| s.as_str());

    let q_sizes_map = q_sizes_file
        .map(pgr::libs::io::read_sizes::<i32>)
        .transpose()?;
    let t_sizes_map = t_sizes_file
        .map(pgr::libs::io::read_sizes::<i32>)
        .transpose()?;

    for line in reader.lines() {
        let line = line?;
        if (line.trim().is_empty() || line.starts_with('#')) && Psl::from_str(&line).is_err() {
            writer.write_fmt(format_args!("{}\n", line))?;
            continue;
        }

        let mut psl: Psl = match line.parse() {
            Ok(p) => p,
            Err(_) => {
                if !line.starts_with('#') && !line.trim().is_empty() {
                    log::warn!("failed to parse psl line, passing through: {}", line);
                }
                writer.write_fmt(format_args!("{}\n", line))?;
                continue;
            }
        };

        if let Some(sizes_map) = &q_sizes_map {
            if !psl.lift_query(sizes_map) {
                log::warn!("failed to lift query: {}", psl.q_name);
            }
        }

        if let Some(sizes_map) = &t_sizes_map {
            if !psl.lift_target(sizes_map) {
                log::warn!("failed to lift target: {}", psl.t_name);
            }
        }

        psl.write_to(&mut writer)?;
    }

    Ok(())
}
