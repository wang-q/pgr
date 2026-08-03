//! Subcommands for runlist interval operations (migrated from `spanr`).

mod compare;
mod cover;
mod coverage;
mod merge;
mod span;

use clap::{ArgMatches, Command};

/// Build the `pgr runlist` subcommand tree.
pub fn make_subcommand() -> Command {
    Command::new("runlist")
        .about("Runlist interval operations (cover, coverage, span, compare, merge)")
        .subcommand(compare::make_subcommand())
        .subcommand(cover::make_subcommand())
        .subcommand(coverage::make_subcommand())
        .subcommand(merge::make_subcommand())
        .subcommand(span::make_subcommand())
}

/// Dispatch `pgr runlist` subcommands.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("compare", sub_matches)) => compare::execute(sub_matches),
        Some(("cover", sub_matches)) => cover::execute(sub_matches),
        Some(("coverage", sub_matches)) => coverage::execute(sub_matches),
        Some(("merge", sub_matches)) => merge::execute(sub_matches),
        Some(("span", sub_matches)) => span::execute(sub_matches),
        _ => unreachable!("runlist subcommand match"),
    }
}
