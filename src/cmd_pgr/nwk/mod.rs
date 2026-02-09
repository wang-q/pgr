use clap::*;

pub mod comment;
pub mod distance;
pub mod indent;
pub mod label;
pub mod order;
pub mod prune;
pub mod rename;
pub mod replace;
pub mod reroot;
pub mod stat;
pub mod subtree;
pub mod to_dot;
pub mod to_forest;
pub mod to_tex;
pub mod topo;
pub mod utils;

pub fn make_subcommand() -> Command {
    Command::new("nwk")
        .about("Manipulate, analyze, and visualize Newick trees")
        .after_help(
            r###"
This suite of tools provides a comprehensive set of operations for phylogenetic trees in Newick format.

Subcommand groups:
* info: stat / label / distance
* ops:  order / prune / rename / replace / reroot / subtree / topo
* viz:  comment / indent

"###,
        )
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(comment::make_subcommand())
        .subcommand(distance::make_subcommand())
        .subcommand(indent::make_subcommand())
        .subcommand(label::make_subcommand())
        .subcommand(order::make_subcommand())
        .subcommand(prune::make_subcommand())
        .subcommand(rename::make_subcommand())
        .subcommand(replace::make_subcommand())
        .subcommand(reroot::make_subcommand())
        .subcommand(stat::make_subcommand())
        .subcommand(subtree::make_subcommand())
        .subcommand(to_dot::make_subcommand())
        .subcommand(to_forest::make_subcommand())
        .subcommand(to_tex::make_subcommand())
        .subcommand(topo::make_subcommand())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("comment", sub_matches)) => comment::execute(sub_matches),
        Some(("distance", sub_matches)) => distance::execute(sub_matches),
        Some(("indent", sub_matches)) => indent::execute(sub_matches),
        Some(("label", sub_matches)) => label::execute(sub_matches),
        Some(("order", sub_matches)) => order::execute(sub_matches),
        Some(("prune", sub_matches)) => prune::execute(sub_matches),
        Some(("rename", sub_matches)) => rename::execute(sub_matches),
        Some(("replace", sub_matches)) => replace::execute(sub_matches),
        Some(("reroot", sub_matches)) => reroot::execute(sub_matches),
        Some(("stat", sub_matches)) => stat::execute(sub_matches),
        Some(("subtree", sub_matches)) => subtree::execute(sub_matches),
        Some(("to-dot", sub_matches)) => to_dot::execute(sub_matches),
        Some(("to-forest", sub_matches)) => to_forest::execute(sub_matches),
        Some(("to-tex", sub_matches)) => to_tex::execute(sub_matches),
        Some(("topo", sub_matches)) => topo::execute(sub_matches),
        _ => unreachable!(),
    }
}
