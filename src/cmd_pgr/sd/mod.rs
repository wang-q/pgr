//! Subcommands for segmental duplication (SD) detection and analysis.

mod align;
mod cluster;
mod cover;
mod decompose;
mod run;
mod search;

use clap::{ArgMatches, Command};

/// Build the `pgr sd` subcommand tree.
pub fn make_subcommand() -> Command {
    Command::new("sd")
        .about("Segmental duplication detection and analysis")
        .subcommand(align::make_subcommand())
        .subcommand(cluster::make_subcommand())
        .subcommand(cover::make_subcommand())
        .subcommand(decompose::make_subcommand())
        .subcommand(run::make_subcommand())
        .subcommand(search::make_subcommand())
}

/// Dispatch `pgr sd` subcommands.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("align", sub_matches)) => align::execute(sub_matches),
        Some(("cluster", sub_matches)) => cluster::execute(sub_matches),
        Some(("cover", sub_matches)) => cover::execute(sub_matches),
        Some(("decompose", sub_matches)) => decompose::execute(sub_matches),
        Some(("run", sub_matches)) => run::execute(sub_matches),
        Some(("search", sub_matches)) => search::execute(sub_matches),
        _ => unreachable!("sd subcommand match"),
    }
}
