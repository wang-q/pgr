use anyhow::Context;
use clap::{ArgMatches, Command};
/// Build the clap subcommand for rc.
pub fn make_subcommand() -> Command {
    Command::new("rc")
        .about("Reverse-complements psl")
        .after_help(
            r###"
Reverse-complement psl.

Examples:
   pgr psl rc in.psl -o out.psl
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg().help("Input PSL file. [stdin] for standard input"))
        .arg(crate::cmd_pgr::args::outfile_arg())
}
/// Execute the rc command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = crate::cmd_pgr::args::get_infile(args);
    let output = crate::cmd_pgr::args::get_outfile(args);

    let reader =
        pgr::reader(input).with_context(|| format!("Failed to open reader for {}", input))?;
    let mut writer =
        pgr::writer(output).with_context(|| format!("Failed to open writer for {}", output))?;

    pgr::libs::fmt::psl::rc_records(reader, &mut writer)?;

    Ok(())
}
