use clap::*;

pub mod indent;
pub mod label;
pub mod stat;
pub mod utils;

pub fn make_subcommand() -> Command {
    Command::new("nwk")
        .about("Newick tools")
        .subcommand(indent::make_subcommand())
        .subcommand(label::make_subcommand())
        .subcommand(stat::make_subcommand())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("indent", sub_matches)) => indent::execute(sub_matches),
        Some(("label", sub_matches)) => label::execute(sub_matches),
        Some(("stat", sub_matches)) => stat::execute(sub_matches),
        _ => unreachable!(),
    }
}
