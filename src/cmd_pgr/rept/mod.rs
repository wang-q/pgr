//! Subcommands for repeat detection and masking.

mod e_align;
mod e_kmer;
mod masker;
mod s_align;
mod s_kmer;
mod trf;

use clap::{ArgMatches, Command};

/// Build the `pgr rept` subcommand tree.
pub fn make_subcommand() -> Command {
    Command::new("rept")
        .about("Detects repetitive regions in a genome")
        .subcommand_required(true)
        .subcommand(e_kmer::make_subcommand())
        .subcommand(e_align::make_subcommand())
        .subcommand(masker::make_subcommand())
        .subcommand(s_align::make_subcommand())
        .subcommand(s_kmer::make_subcommand())
        .subcommand(trf::make_subcommand())
}

/// Dispatch `pgr rept` subcommands.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("e-kmer", sub_matches)) => e_kmer::execute(sub_matches),
        Some(("e-align", sub_matches)) => e_align::execute(sub_matches),
        Some(("masker", sub_matches)) => masker::execute(sub_matches),
        Some(("s-align", sub_matches)) => s_align::execute(sub_matches),
        Some(("s-kmer", sub_matches)) => s_kmer::execute(sub_matches),
        Some(("trf", sub_matches)) => trf::execute(sub_matches),
        _ => unreachable!("rept subcommand match"),
    }
}
