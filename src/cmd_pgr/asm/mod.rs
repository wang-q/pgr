pub mod cns;
pub mod common;
pub mod contig;
pub mod layout;
pub mod map;
pub mod olc;
pub mod ovlp;
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
        .subcommand(ovlp::make_subcommand())
        .subcommand(layout::make_subcommand())
        .subcommand(cns::make_subcommand())
        .subcommand(olc::make_subcommand())
        .subcommand(map::make_subcommand())
}
/// Execute the asm command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("contig", sub_matches)) => contig::execute(sub_matches),
        Some(("unitig", sub_matches)) => unitig::execute(sub_matches),
        Some(("ovlp", sub_matches)) => ovlp::execute(sub_matches),
        Some(("layout", sub_matches)) => layout::execute(sub_matches),
        Some(("cns", sub_matches)) => cns::execute(sub_matches),
        Some(("olc", sub_matches)) => olc::execute(sub_matches),
        Some(("map", sub_matches)) => map::execute(sub_matches),
        _ => Ok(()),
    }
}
