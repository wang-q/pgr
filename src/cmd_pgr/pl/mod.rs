pub mod p2m;
pub mod trf;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("pl")
        .about("Pipeline tools")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(p2m::make_subcommand())
        .subcommand(trf::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("p2m", sub_matches)) => p2m::execute(sub_matches),
        Some(("trf", sub_matches)) => trf::execute(sub_matches),
        _ => unreachable!(),
    }
}
