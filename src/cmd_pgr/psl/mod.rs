pub mod histo;
pub mod stats;
pub mod swap;
pub mod tochain;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("psl")
        .about("Psl tools")
        .subcommand(histo::make_subcommand())
        .subcommand(stats::make_subcommand())
        .subcommand(swap::make_subcommand())
        .subcommand(tochain::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("histo", sub_matches)) => histo::execute(sub_matches),
        Some(("stats", sub_matches)) => stats::execute(sub_matches),
        Some(("swap", sub_matches)) => swap::execute(sub_matches),
        Some(("tochain", sub_matches)) => tochain::execute(sub_matches),
        _ => Ok(()),
    }
}
