pub mod split;
pub mod stitch;

use clap::{Command, ArgMatches};

pub fn make_subcommand() -> Command {
    Command::new("chain")
        .about("Chain tools")
        .subcommand_required(true)
        .subcommand(split::make_subcommand())
        .subcommand(stitch::make_subcommand())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("split", sub_matches)) => split::execute(sub_matches),
        Some(("stitch", sub_matches)) => stitch::execute(sub_matches),
        _ => unreachable!(),
    }
}
