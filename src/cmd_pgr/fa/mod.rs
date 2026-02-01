pub mod to2bit;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("fa")
        .about("Fasta tools")
        .subcommand(to2bit::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("to2bit", sub_matches)) => to2bit::execute(sub_matches),
        _ => Ok(()),
    }
}
