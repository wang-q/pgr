pub fn make_subcommand() -> clap::Command {
    clap::Command::new("chaining")
        .about("Chaining alignment blocks")
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        _ => Ok(()),
    }
}
