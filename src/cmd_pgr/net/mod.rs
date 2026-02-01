pub mod syntenic;
pub mod toaxt;
pub mod split;
pub mod subset;

use clap::{ArgMatches, Command};

pub fn make_subcommand() -> Command {
    Command::new("net")
        .about("Net tools")
        .subcommand_required(true)
        .subcommand(syntenic::make_subcommand())
        .subcommand(toaxt::make_subcommand())
        .subcommand(split::make_subcommand())
        .subcommand(subset::make_subcommand())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("syntenic", sub_matches)) => syntenic::execute(sub_matches),
        Some(("to-axt", sub_matches)) => toaxt::execute(sub_matches),
        Some(("split", sub_matches)) => split::execute(sub_matches),
        Some(("subset", sub_matches)) => subset::execute(sub_matches),
        _ => unreachable!(),
    }
}
