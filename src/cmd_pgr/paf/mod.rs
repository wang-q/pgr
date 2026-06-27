pub mod graph;
pub mod index;
pub mod query;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("paf")
        .about("Manipulate PAF (Pairwise mApping Format) files")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(index::make_subcommand())
        .subcommand(query::make_subcommand())
        .subcommand(graph::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("index", sub_matches)) => index::execute(sub_matches),
        Some(("query", sub_matches)) => query::execute(sub_matches),
        Some(("graph", sub_matches)) => graph::execute(sub_matches),
        _ => Ok(()),
    }
}
