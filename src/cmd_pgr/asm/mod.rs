pub mod contig;
pub mod map;
pub mod unitig;

use clap::{ArgMatches, Command};
/// Build the clap subcommand for asm.
pub fn make_subcommand() -> Command {
    Command::new("asm")
        .about("Assembles reads into contigs/unitigs and maps reads back")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(contig::make_subcommand())
        .subcommand(unitig::make_subcommand())
        .subcommand(map::make_subcommand())
}
/// Execute the asm command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("contig", sub_matches)) => contig::execute(sub_matches),
        Some(("unitig", sub_matches)) => unitig::execute(sub_matches),
        Some(("map", sub_matches)) => map::execute(sub_matches),
        _ => Ok(()),
    }
}
