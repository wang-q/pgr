pub mod histo;
pub mod tochain;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("psl")
        .about("Psl tools")
        .subcommand(histo::make_subcommand())
        .subcommand(tochain::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("histo", sub_matches)) => histo::execute(sub_matches),
        Some(("tochain", sub_matches)) => tochain::execute(sub_matches),
        _ => Ok(()),
    }
}
