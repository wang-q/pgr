extern crate clap;
use clap::*;

mod cmd_pgr;

fn main() -> anyhow::Result<()> {
    let app = Command::new("pgr")
        .version(crate_version!())
        .author(crate_authors!())
        .about("`pgr` - Practical Genome Refiner")
        .propagate_version(true)
        .arg_required_else_help(true)
        .color(ColorChoice::Auto)
        .subcommand(cmd_pgr::pipeline::make_subcommand())
        .subcommand(cmd_pgr::ir::make_subcommand())
        .subcommand(cmd_pgr::rept::make_subcommand())
        .subcommand(cmd_pgr::trf::make_subcommand())
        .subcommand(cmd_pgr::ms2dna::make_subcommand())
        .subcommand(cmd_pgr::axt::make_subcommand())
        .subcommand(cmd_pgr::chain::make_subcommand())
        .subcommand(cmd_pgr::lav::make_subcommand())
        .subcommand(cmd_pgr::net::make_subcommand())
        .subcommand(cmd_pgr::psl::make_subcommand())
        .subcommand(cmd_pgr::twobit::make_subcommand())
        .subcommand(cmd_pgr::fa::make_subcommand())
        .subcommand(cmd_pgr::fq::make_subcommand())
        .after_help(
            r###"Subcommand groups:

* Fasta files
    * info: size / count / masked / n50
    * records: one / some / order / split
    * transform: replace / rc / filter / dedup / mask / sixframe
    * indexing: gz / range / prefilter

* Genome alignments:
    * chain
    * net
    * axt
    * lav
    * psl
    * 2bit
    * fa
    * fq

* Repeats:
    * ir / rept / trf

"###,
        );

    // Check which subcomamnd the user ran...
    match app.get_matches().subcommand() {
        Some(("pipeline", sub_matches)) => cmd_pgr::pipeline::execute(sub_matches),
        Some(("ir", sub_matches)) => cmd_pgr::ir::execute(sub_matches),
        Some(("rept", sub_matches)) => cmd_pgr::rept::execute(sub_matches),
        Some(("trf", sub_matches)) => cmd_pgr::trf::execute(sub_matches),
        Some(("ms2dna", sub_matches)) => cmd_pgr::ms2dna::execute(sub_matches),
        Some(("axt", sub_matches)) => cmd_pgr::axt::execute(sub_matches),
        Some(("chain", sub_matches)) => cmd_pgr::chain::execute(sub_matches),
        Some(("lav", sub_matches)) => cmd_pgr::lav::execute(sub_matches),
        Some(("net", sub_matches)) => cmd_pgr::net::execute(sub_matches),
        Some(("psl", sub_matches)) => cmd_pgr::psl::execute(sub_matches),
        Some(("2bit", sub_matches)) => cmd_pgr::twobit::execute(sub_matches),
        Some(("fa", sub_matches)) => cmd_pgr::fa::execute(sub_matches),
        Some(("fq", sub_matches)) => cmd_pgr::fq::execute(sub_matches),
        _ => unreachable!(),
    }?;

    Ok(())
}

// TODO: paralog
