use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};
use std::io::Write;

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
* Lines that fail to parse as PSL records are skipped with a warning.
  Use --strict to turn parse failures into hard errors.

Examples:
1. Lift query coordinates:
   pgr psl lift input.psl --q-sizes chrom.sizes > output.psl

2. Lift target coordinates:
   pgr psl lift input.psl --t-sizes chrom.sizes > output.psl

3. Lift both:
   pgr psl lift input.psl --q-sizes q.sizes --t-sizes t.sizes > output.psl

4. Strict mode (fail on parse errors):
   pgr psl lift input.psl --q-sizes q.sizes --strict -o output.psl
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input PSL file. [stdin] for standard input",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(crate::cmd_pgr::args::q_sizes_arg())
        .arg(crate::cmd_pgr::args::t_sizes_arg())
        .arg(
            Arg::new("strict")
                .long("strict")
                .action(ArgAction::SetTrue)
                .help("Fail on parse errors instead of skipping malformed lines"),
        )
}

/// Execute the lift command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let output = crate::cmd_pgr::args::get_outfile(args);
    let strict = args.get_flag("strict");

    let reader =
        pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;
    let mut writer =
        pgr::writer(output).with_context(|| format!("Failed to open writer for {}", output))?;

    let q_sizes_file = args.get_one::<String>("q_sizes").map(|s| s.as_str());
    let t_sizes_file = args.get_one::<String>("t_sizes").map(|s| s.as_str());

    let q_sizes_map = q_sizes_file
        .map(pgr::libs::io::read_sizes::<i32>)
        .transpose()?;
    let t_sizes_map = t_sizes_file
        .map(pgr::libs::io::read_sizes::<i32>)
        .transpose()?;

    pgr::libs::fmt::psl::lift(
        reader,
        &mut writer,
        q_sizes_map.as_ref(),
        t_sizes_map.as_ref(),
        strict,
    )?;

    writer.flush()?;
    Ok(())
}
