pub mod lastz;
pub mod to_psl;

use clap::{ArgMatches, Command};

pub fn make_subcommand() -> Command {
    Command::new("lav")
        .about("Manipulate LAV alignment files")
        .subcommand_required(true)
        .subcommand(lastz::make_subcommand())
        .subcommand(to_psl::make_subcommand())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("lastz", sub_matches)) => lastz::execute(sub_matches),
        Some(("to-psl", sub_matches)) => to_psl::execute(sub_matches),
        _ => unreachable!(),
    }
}
