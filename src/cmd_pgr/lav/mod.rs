pub mod to_psl;

use clap::{ArgMatches, Command};

pub fn make_subcommand() -> Command {
    Command::new("lav")
        .about("LAV tools")
        .subcommand_required(true)
        .subcommand(to_psl::make_subcommand())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("to-psl", sub_matches)) => to_psl::execute(sub_matches),
        _ => unreachable!(),
    }
}
