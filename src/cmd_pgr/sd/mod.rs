//! Subcommands for segmental duplication (SD) detection and analysis.

mod search;

use clap::{ArgMatches, Command};

/// Build the `pgr sd` subcommand tree.
pub fn make_subcommand() -> Command {
    Command::new("sd")
        .about("Segmental duplication detection and analysis")
        .subcommand(search::make_subcommand())
}

/// Dispatch `pgr sd` subcommands.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("search", sub_matches)) => search::execute(sub_matches),
        _ => unreachable!("sd subcommand match"),
    }
}
