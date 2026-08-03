//! Subcommands for runlist interval operations (migrated from `spanr`).

mod combine;
mod compare;
mod convert;
mod cover;
mod coverage;
mod genome;
mod merge;
mod some;
mod span;
mod split;
mod stat;
mod statop;

use clap::{ArgMatches, Command};

/// Build the `pgr runlist` subcommand tree.
pub fn make_subcommand() -> Command {
    Command::new("runlist")
        .about("Runlist interval operations")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(combine::make_subcommand())
        .subcommand(compare::make_subcommand())
        .subcommand(convert::make_subcommand())
        .subcommand(cover::make_subcommand())
        .subcommand(coverage::make_subcommand())
        .subcommand(genome::make_subcommand())
        .subcommand(merge::make_subcommand())
        .subcommand(some::make_subcommand())
        .subcommand(span::make_subcommand())
        .subcommand(split::make_subcommand())
        .subcommand(stat::make_subcommand())
        .subcommand(statop::make_subcommand())
}

/// Dispatch `pgr runlist` subcommands.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("combine", sub_matches)) => combine::execute(sub_matches),
        Some(("compare", sub_matches)) => compare::execute(sub_matches),
        Some(("convert", sub_matches)) => convert::execute(sub_matches),
        Some(("cover", sub_matches)) => cover::execute(sub_matches),
        Some(("coverage", sub_matches)) => coverage::execute(sub_matches),
        Some(("genome", sub_matches)) => genome::execute(sub_matches),
        Some(("merge", sub_matches)) => merge::execute(sub_matches),
        Some(("some", sub_matches)) => some::execute(sub_matches),
        Some(("span", sub_matches)) => span::execute(sub_matches),
        Some(("split", sub_matches)) => split::execute(sub_matches),
        Some(("stat", sub_matches)) => stat::execute(sub_matches),
        Some(("statop", sub_matches)) => statop::execute(sub_matches),
        _ => unreachable!("runlist subcommand match"),
    }
}
