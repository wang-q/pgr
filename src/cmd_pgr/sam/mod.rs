pub mod to_rg;

use clap::{ArgMatches, Command};

/// Build the clap subcommand for sam.
pub fn make_subcommand() -> Command {
    Command::new("sam")
        .about("Manipulates SAM alignment files")
        .subcommand(to_rg::make_subcommand())
}

/// Execute the sam command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("to-rg", sub_matches)) => to_rg::execute(sub_matches),
        _ => Ok(()),
    }
}
