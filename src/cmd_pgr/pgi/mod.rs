//! Subcommands for pgr genome index (.pgi) files.

mod build;
mod stat;
mod to_hv;

use clap::{ArgMatches, Command};

/// Build the `pgr pgi` subcommand tree.
pub fn make_subcommand() -> Command {
    Command::new("pgi")
        .about("Manages pgr genome index (.pgi) files")
        .subcommand_required(true)
        .subcommand(build::make_subcommand())
        .subcommand(stat::make_subcommand())
        .subcommand(to_hv::make_subcommand())
}

/// Dispatch `pgr pgi` subcommands.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("build", sub_matches)) => build::execute(sub_matches),
        Some(("stat", sub_matches)) => stat::execute(sub_matches),
        Some(("to-hv", sub_matches)) => to_hv::execute(sub_matches),
        _ => unreachable!("pgi subcommand match"),
    }
}
