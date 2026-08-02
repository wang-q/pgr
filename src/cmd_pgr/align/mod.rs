//! Subcommands for pairwise genome alignment.

mod pgi;

use clap::{ArgMatches, Command};

/// Build the `pgr align` subcommand tree.
pub fn make_subcommand() -> Command {
    Command::new("align")
        .about("Aligns genomes or .pgi indexes into PSL blocks")
        .subcommand(pgi::make_subcommand())
}

/// Dispatch `pgr align` subcommands.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("pgi", sub_matches)) => pgi::execute(sub_matches),
        _ => unreachable!("align subcommand match"),
    }
}
