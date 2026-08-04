//! Subcommands for `pgr rg` — line-oriented operations on `.rg` range files.

mod count;
mod cover;
mod coverage;
mod prop;
mod sort;

use clap::{ArgMatches, Command};

/// Build the `pgr rg` subcommand tree.
pub fn make_subcommand() -> Command {
    Command::new("rg")
        .about("Operates on .rg range lines")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(cover::make_subcommand())
        .subcommand(count::make_subcommand())
        .subcommand(coverage::make_subcommand())
        .subcommand(prop::make_subcommand())
        .subcommand(sort::make_subcommand())
}

/// Dispatch `pgr rg` subcommands.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("cover", sub_matches)) => cover::execute(sub_matches),
        Some(("count", sub_matches)) => count::execute(sub_matches),
        Some(("coverage", sub_matches)) => coverage::execute(sub_matches),
        Some(("prop", sub_matches)) => prop::execute(sub_matches),
        Some(("sort", sub_matches)) => sort::execute(sub_matches),
        _ => unreachable!("rg subcommand match"),
    }
}
