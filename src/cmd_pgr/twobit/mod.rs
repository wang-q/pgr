pub mod masked;
pub mod range;
pub mod size;
pub mod some;
pub mod to_fa;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("2bit")
        .about("Manage 2bit files")
        .after_help(
            r###"Subcommand groups:

* info: masked / size
* subset: range / some
* transform: to-fa

"###,
        )
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(masked::make_subcommand())
        .subcommand(range::make_subcommand())
        .subcommand(size::make_subcommand())
        .subcommand(some::make_subcommand())
        .subcommand(to_fa::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("masked", sub_matches)) => masked::execute(sub_matches),
        Some(("range", sub_matches)) => range::execute(sub_matches),
        Some(("size", sub_matches)) => size::execute(sub_matches),
        Some(("some", sub_matches)) => some::execute(sub_matches),
        Some(("to-fa", sub_matches)) => to_fa::execute(sub_matches),
        _ => Ok(()),
    }
}
