pub mod cover;
pub mod name;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("fas")
        .about("Block FA tools")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(cover::make_subcommand())
        .subcommand(name::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("cover", sub_matches)) => cover::execute(sub_matches),
        Some(("name", sub_matches)) => name::execute(sub_matches),
        _ => unreachable!(),
    }
}
