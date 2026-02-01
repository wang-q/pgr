pub mod class;
pub mod filter;
pub mod split;
pub mod subset;
pub mod syntenic;
pub mod toaxt;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("net")
        .about("Net tools")
        .subcommand(class::make_subcommand())
        .subcommand(filter::make_subcommand())
        .subcommand(split::make_subcommand())
        .subcommand(subset::make_subcommand())
        .subcommand(syntenic::make_subcommand())
        .subcommand(toaxt::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("class", sub_matches)) => class::execute(sub_matches),
        Some(("filter", sub_matches)) => filter::execute(sub_matches),
        Some(("split", sub_matches)) => split::execute(sub_matches),
        Some(("subset", sub_matches)) => subset::execute(sub_matches),
        Some(("syntenic", sub_matches)) => syntenic::execute(sub_matches),
        Some(("to-axt", sub_matches)) => toaxt::execute(sub_matches),
        _ => Ok(()),
    }
}
