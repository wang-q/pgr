use anyhow::Context;
use clap::{Arg, ArgMatches, Command};
use std::io::Write;
/// Build the clap subcommand for swap.
pub fn make_subcommand() -> Command {
    Command::new("swap")
        .about("Reverses target and query in psls")
        .after_help(
            r###"
Reverse target and query in psls.

Examples:
   pgr psl swap in.psl -o out.psl
   pgr psl swap --no-rc in.psl -o out.psl
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg().help("Input PSL file. [stdin] for standard input"))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("no_rc")
                .short('n')
                .long("no-rc")
                .action(clap::ArgAction::SetTrue)
                .help("Don't reverse-complement PSL if needed, instead make target strand explict"),
        )
}
/// Execute the swap command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = crate::cmd_pgr::args::get_infile(args);
    let output = crate::cmd_pgr::args::get_outfile(args);
    let no_rc = args.get_flag("no_rc");

    let reader =
        pgr::reader(input).with_context(|| format!("Failed to open reader for {}", input))?;
    let mut writer =
        pgr::writer(output).with_context(|| format!("Failed to open writer for {}", output))?;

    pgr::libs::fmt::psl::swap_records(reader, &mut writer, no_rc)?;

    writer.flush()?;
    Ok(())
}
