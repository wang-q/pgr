pub mod common;
pub mod frac;
pub mod hv;
pub mod mash;
pub mod mini;
pub mod pgi;

use clap::{ArgMatches, Command};
/// Build the clap subcommand for dist.
pub fn make_subcommand() -> Command {
    Command::new("dist")
        .about("Computes distance/similarity metrics")
        .after_help(
            r###"Subcommand groups:

* sketch distances: mini (minimizer, ranking) / mash (MinHash, Mash-compatible)
  / frac (FracMinHash, unbiased numeric ANI with CI)
* other: hv (hypervectors) / pgi (syncmer index merge)

"###,
        )
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(hv::make_subcommand())
        .subcommand(pgi::make_subcommand())
        .subcommand(mini::make_subcommand())
        .subcommand(mash::make_subcommand())
        .subcommand(frac::make_subcommand())
}
/// Execute the dist command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("hv", sub_matches)) => hv::execute(sub_matches),
        Some(("pgi", sub_matches)) => pgi::execute(sub_matches),
        Some(("mini", sub_matches)) => mini::execute(sub_matches),
        Some(("mash", sub_matches)) => mash::execute(sub_matches),
        Some(("frac", sub_matches)) => frac::execute(sub_matches),
        _ => Ok(()),
    }
}
