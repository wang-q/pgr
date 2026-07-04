pub mod cc;
pub mod cut;
pub mod dbscan;
pub mod eval;
pub mod hier;
pub mod k_medoids;
pub mod mcl;
pub mod nj;
pub mod upgma;

use clap::{ArgMatches, Command};
/// Build the clap subcommand for clust.
pub fn make_subcommand() -> Command {
    Command::new("clust")
        .about("Clusters entries via various algorithms")
        .after_help(
            r###"Subcommand groups:

* Tree: hier, nj, upgma
* Flat: cc, cut, dbscan, k-medoids, mcl
* Eval: eval

"###,
        )
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(cc::make_subcommand())
        .subcommand(cut::make_subcommand())
        .subcommand(dbscan::make_subcommand())
        .subcommand(eval::make_subcommand())
        .subcommand(hier::make_subcommand())
        .subcommand(k_medoids::make_subcommand())
        .subcommand(mcl::make_subcommand())
        .subcommand(nj::make_subcommand())
        .subcommand(upgma::make_subcommand())
}
/// Execute the clust command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("cc", sub_matches)) => cc::execute(sub_matches),
        Some(("cut", sub_matches)) => cut::execute(sub_matches),
        Some(("dbscan", sub_matches)) => dbscan::execute(sub_matches),
        Some(("eval", sub_matches)) => eval::execute(sub_matches),
        Some(("hier", sub_matches)) => hier::execute(sub_matches),
        Some(("k-medoids", sub_matches)) => k_medoids::execute(sub_matches),
        Some(("mcl", sub_matches)) => mcl::execute(sub_matches),
        Some(("nj", sub_matches)) => nj::execute(sub_matches),
        Some(("upgma", sub_matches)) => upgma::execute(sub_matches),
        _ => Ok(()),
    }
}
