use clap::*;

pub mod distance;
pub mod indent;
pub mod label;
pub mod rename;
pub mod replace;
pub mod stat;
pub mod topo;
pub mod utils;

pub fn make_subcommand() -> Command {
    Command::new("nwk")
        .about("Newick tools")
        .after_help(
            r###"Subcommand groups:

* info: stat / label / distance
* ops:  rename / replace / topo
* viz:  indent

"###,
        )
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(distance::make_subcommand())
        .subcommand(indent::make_subcommand())
        .subcommand(label::make_subcommand())
        .subcommand(rename::make_subcommand())
        .subcommand(replace::make_subcommand())
        .subcommand(stat::make_subcommand())
        .subcommand(topo::make_subcommand())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("distance", sub_matches)) => distance::execute(sub_matches),
        Some(("indent", sub_matches)) => indent::execute(sub_matches),
        Some(("label", sub_matches)) => label::execute(sub_matches),
        Some(("rename", sub_matches)) => rename::execute(sub_matches),
        Some(("replace", sub_matches)) => replace::execute(sub_matches),
        Some(("stat", sub_matches)) => stat::execute(sub_matches),
        Some(("topo", sub_matches)) => topo::execute(sub_matches),
        _ => unreachable!(),
    }
}
