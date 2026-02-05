pub mod to_dna;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("ms")
        .about("Hudson's ms simulator tools")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(to_dna::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("to-dna", sub_matches)) => to_dna::execute(sub_matches),
        _ => Ok(()),
    }
}
