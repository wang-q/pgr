pub mod hv;
pub mod seq;
pub mod vector;

use clap::{ArgMatches, Command};
/// Build the clap subcommand for dist.
pub fn make_subcommand() -> Command {
    Command::new("dist")
        .about("Distance/Similarity metrics")
        .after_help(
            r###"Subcommand groups:

* distance: hv / seq / vector

"###,
        )
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(hv::make_subcommand())
        .subcommand(seq::make_subcommand())
        .subcommand(vector::make_subcommand())
}
/// Execute the dist command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("hv", sub_matches)) => hv::execute(sub_matches),
        Some(("seq", sub_matches)) => seq::execute(sub_matches),
        Some(("vector", sub_matches)) => vector::execute(sub_matches),
        _ => Ok(()),
    }
}
