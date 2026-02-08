use clap::*;

pub mod stat;

pub fn make_subcommand() -> Command {
    Command::new("nwk")
        .about("Newick tools")
        .subcommand(stat::make_subcommand())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("stat", sub_matches)) => stat::execute(sub_matches),
        _ => unreachable!(),
    }
}
