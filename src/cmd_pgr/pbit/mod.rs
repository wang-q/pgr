pub mod append;
pub mod create;
pub mod range;
pub mod some;
pub mod stat;
pub mod to_fa;

use clap::{ArgMatches, Command};

/// Build the clap subcommand for pbit.
pub fn make_subcommand() -> Command {
    Command::new("pbit")
        .about("Manages pbit (population 2bit + delta) files")
        .after_help(
            r###"Subcommand groups:

* build:     create / append
* info:      stat
* subset:    range / some
* transform: to-fa

"###,
        )
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(create::make_subcommand())
        .subcommand(append::make_subcommand())
        .subcommand(range::make_subcommand())
        .subcommand(some::make_subcommand())
        .subcommand(stat::make_subcommand())
        .subcommand(to_fa::make_subcommand())
}

/// Execute the pbit command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("create", sub_matches)) => create::execute(sub_matches),
        Some(("append", sub_matches)) => append::execute(sub_matches),
        Some(("range", sub_matches)) => range::execute(sub_matches),
        Some(("some", sub_matches)) => some::execute(sub_matches),
        Some(("stat", sub_matches)) => stat::execute(sub_matches),
        Some(("to-fa", sub_matches)) => to_fa::execute(sub_matches),
        _ => Ok(()),
    }
}
