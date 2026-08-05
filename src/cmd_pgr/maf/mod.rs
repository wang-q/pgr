pub mod to_fas;
pub mod to_paf;

use clap::{ArgMatches, Command};
/// Build the clap subcommand for maf.
pub fn make_subcommand() -> Command {
    Command::new("maf")
        .about("Manipulates MAF alignment files")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(to_fas::make_subcommand())
        .subcommand(to_paf::make_subcommand())
}
/// Execute the maf command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("to-fas", sub_matches)) => to_fas::execute(sub_matches),
        Some(("to-paf", sub_matches)) => to_paf::execute(sub_matches),
        _ => Ok(()),
    }
}
