pub mod to_psl;

use clap::{ArgMatches, Command};
/// Build the clap subcommand for lav.
pub fn make_subcommand() -> Command {
    Command::new("lav")
        .about("Manipulates LAV alignment files")
        .subcommand_required(true)
        .subcommand(to_psl::make_subcommand())
}
/// Execute the lav command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("to-psl", sub_matches)) => to_psl::execute(sub_matches),
        _ => Ok(()),
    }
}
