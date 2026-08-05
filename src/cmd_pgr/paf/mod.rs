pub mod graph;
pub mod index;
pub mod query;
pub mod stat;
pub mod to_bed;
pub mod to_fas;
pub mod to_gfa;
pub mod to_maf;
pub mod to_vcf;

use clap::{ArgMatches, Command};
/// Build the clap subcommand for paf.
pub fn make_subcommand() -> Command {
    Command::new("paf")
        .about("Manipulates PAF (Pairwise mApping Format) files")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(index::make_subcommand())
        .subcommand(query::make_subcommand())
        .subcommand(to_bed::make_subcommand())
        .subcommand(to_fas::make_subcommand())
        .subcommand(to_maf::make_subcommand())
        .subcommand(to_vcf::make_subcommand())
        .subcommand(to_gfa::make_subcommand())
        .subcommand(graph::make_subcommand())
        .subcommand(stat::make_subcommand())
}
/// Execute the paf command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("index", sub_matches)) => index::execute(sub_matches),
        Some(("query", sub_matches)) => query::execute(sub_matches),
        Some(("to-bed", sub_matches)) => to_bed::execute(sub_matches),
        Some(("to-fas", sub_matches)) => to_fas::execute(sub_matches),
        Some(("to-maf", sub_matches)) => to_maf::execute(sub_matches),
        Some(("to-vcf", sub_matches)) => to_vcf::execute(sub_matches),
        Some(("to-gfa", sub_matches)) => to_gfa::execute(sub_matches),
        Some(("graph", sub_matches)) => graph::execute(sub_matches),
        Some(("stat", sub_matches)) => stat::execute(sub_matches),
        _ => Ok(()),
    }
}
