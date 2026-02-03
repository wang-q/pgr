pub mod concat;
pub mod cover;
pub mod link;
pub mod name;
pub mod subset;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("fas")
        .about("Block FA tools")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(concat::make_subcommand())
        .subcommand(cover::make_subcommand())
        .subcommand(link::make_subcommand())
        .subcommand(name::make_subcommand())
        .subcommand(subset::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("concat", sub_matches)) => concat::execute(sub_matches),
        Some(("cover", sub_matches)) => cover::execute(sub_matches),
        Some(("link", sub_matches)) => link::execute(sub_matches),
        Some(("name", sub_matches)) => name::execute(sub_matches),
        Some(("subset", sub_matches)) => subset::execute(sub_matches),
        _ => unreachable!(),
    }
}
