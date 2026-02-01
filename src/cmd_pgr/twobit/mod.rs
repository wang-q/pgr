pub mod to2bit;
pub mod tofa;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("2bit")
        .about("2bit tools")
        .subcommand(to2bit::make_subcommand())
        .subcommand(tofa::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("to2bit", sub_matches)) => to2bit::execute(sub_matches),
        Some(("tofa", sub_matches)) => tofa::execute(sub_matches),
        _ => Ok(()),
    }
}
