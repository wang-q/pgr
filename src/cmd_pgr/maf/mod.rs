pub mod to_fas;
pub mod to_paf;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("maf")
        .about("Manipulates MAF alignment files")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(to_fas::make_subcommand())
        .subcommand(to_paf::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("to-fas", sub_matches)) => to_fas::execute(sub_matches),
        Some(("to-paf", sub_matches)) => to_paf::execute(sub_matches),
        _ => Ok(()),
    }
}
