pub mod cc;
pub mod dbscan;
pub mod kmedoids;
pub mod mcl;

pub fn make_subcommand() -> clap::Command {
    clap::Command::new("clust")
        .about("Clustering operations")
        .after_help(
            r###"Subcommand groups:

* clustering: cc, dbscan, k-medoids, mcl

"###,
        )
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(cc::make_subcommand())
        .subcommand(dbscan::make_subcommand())
        .subcommand(kmedoids::make_subcommand())
        .subcommand(mcl::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("cc", sub_matches)) => cc::execute(sub_matches),
        Some(("dbscan", sub_matches)) => dbscan::execute(sub_matches),
        Some(("k-medoids", sub_matches)) => kmedoids::execute(sub_matches),
        Some(("mcl", sub_matches)) => mcl::execute(sub_matches),
        _ => Ok(()),
    }
}
