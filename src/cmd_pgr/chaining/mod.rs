pub mod psl;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("chaining")
        .about("Chaining alignment blocks")
        .subcommand(psl::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("psl", sub_matches)) => psl::execute(sub_matches),
        _ => Ok(()),
    }
}
