// Subcommand modules for the `gff` command.
pub mod rg;

use clap::{ArgMatches, Command};
/// Build the clap subcommand for gff.
pub fn make_subcommand() -> Command {
    Command::new("gff")
        .about("Manipulates GFF files")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(rg::make_subcommand())
}
/// Execute the gff command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("rg", sub_matches)) => rg::execute(sub_matches),
        _ => Ok(()),
    }
}
