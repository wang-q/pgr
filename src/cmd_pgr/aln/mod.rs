pub mod expand;
pub mod stat;
pub mod to_paf;
pub mod to_psl;

use clap::{ArgMatches, Command};

/// Build the clap subcommand for 1aln.
pub fn make_subcommand() -> Command {
    Command::new("1aln")
        .about("Reads FastGA .1aln (ONEcode trace-point) alignment files")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(stat::make_subcommand())
        .subcommand(to_paf::make_subcommand())
        .subcommand(to_psl::make_subcommand())
}

/// Execute the 1aln command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("stat", sub_matches)) => stat::execute(sub_matches),
        Some(("to-paf", sub_matches)) => to_paf::execute(sub_matches),
        Some(("to-psl", sub_matches)) => to_psl::execute(sub_matches),
        _ => Ok(()),
    }
}
