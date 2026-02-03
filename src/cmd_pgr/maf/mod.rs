pub mod tofas;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("maf")
        .about("Maf tools")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(tofas::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("tofas", sub_matches)) => tofas::execute(sub_matches),
        _ => Ok(()),
    }
}
