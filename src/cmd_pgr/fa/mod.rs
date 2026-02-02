pub mod masked;
pub mod n50;
pub mod size;
pub mod some;
pub mod to2bit;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("fa")
        .about("Fasta tools")
        .after_help(
            r###"Notes:
* Supports both plain text and gzipped (.gz) files

"###,
        )
        .subcommand(masked::make_subcommand())
        .subcommand(n50::make_subcommand())
        .subcommand(size::make_subcommand())
        .subcommand(some::make_subcommand())
        .subcommand(to2bit::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("masked", sub_matches)) => masked::execute(sub_matches),
        Some(("n50", sub_matches)) => n50::execute(sub_matches),
        Some(("size", sub_matches)) => size::execute(sub_matches),
        Some(("some", sub_matches)) => some::execute(sub_matches),
        Some(("to2bit", sub_matches)) => to2bit::execute(sub_matches),
        _ => Ok(()),
    }
}
