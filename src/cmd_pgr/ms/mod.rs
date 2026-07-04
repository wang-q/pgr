pub mod to_dna;

use clap::{ArgMatches, Command};
/// Build the clap subcommand for ms.
pub fn make_subcommand() -> Command {
    Command::new("ms")
        .about("Hudson's ms simulator tools")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(to_dna::make_subcommand())
}
/// Execute the ms command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("to-dna", sub_matches)) => to_dna::execute(sub_matches),
        _ => Ok(()),
    }
}
