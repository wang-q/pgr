pub mod clean;
pub mod clump;
pub mod ec_kmer;
pub mod ec_overlap;
pub mod extend;
pub mod filter;
pub mod interleave;
pub mod merge;
pub mod norm;
pub mod range;
pub mod s_filter;
pub mod sample;
pub mod split;
pub mod to_fa;
pub mod trim_qual;

use clap::{ArgMatches, Command};
/// Build the clap subcommand for fq.
pub fn make_subcommand() -> Command {
    Command::new("fq")
        .about("Manipulates FASTQ files")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(to_fa::make_subcommand())
        .subcommand(clump::make_subcommand())
        .subcommand(interleave::make_subcommand())
        .subcommand(merge::make_subcommand())
        .subcommand(norm::make_subcommand())
        .subcommand(range::make_subcommand())
        .subcommand(sample::make_subcommand())
        .subcommand(split::make_subcommand())
        .subcommand(clean::make_subcommand())
        .subcommand(ec_kmer::make_subcommand())
        .subcommand(ec_overlap::make_subcommand())
        .subcommand(extend::make_subcommand())
        .subcommand(filter::make_subcommand())
        .subcommand(s_filter::make_subcommand())
        .subcommand(trim_qual::make_subcommand())
}
/// Execute the fq command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("to-fa", sub_matches)) => to_fa::execute(sub_matches),
        Some(("clump", sub_matches)) => clump::execute(sub_matches),
        Some(("interleave", sub_matches)) | Some(("il", sub_matches)) => {
            interleave::execute(sub_matches)
        }
        Some(("merge", sub_matches)) => merge::execute(sub_matches),
        Some(("norm", sub_matches)) => norm::execute(sub_matches),
        Some(("range", sub_matches)) => range::execute(sub_matches),
        Some(("sample", sub_matches)) => sample::execute(sub_matches),
        Some(("split", sub_matches)) => split::execute(sub_matches),
        Some(("clean", sub_matches)) => clean::execute(sub_matches),
        Some(("ec-kmer", sub_matches)) => ec_kmer::execute(sub_matches),
        Some(("ec-overlap", sub_matches)) => ec_overlap::execute(sub_matches),
        Some(("extend", sub_matches)) => extend::execute(sub_matches),
        Some(("filter", sub_matches)) => filter::execute(sub_matches),
        Some(("s-filter", sub_matches)) => s_filter::execute(sub_matches),
        Some(("trim-qual", sub_matches)) => trim_qual::execute(sub_matches),
        _ => Ok(()),
    }
}
