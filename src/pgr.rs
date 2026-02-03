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
        .subcommand(cmd_pgr::ms2dna::make_subcommand())
        .subcommand(cmd_pgr::axt::make_subcommand())
        .subcommand(cmd_pgr::chain::make_subcommand())
        .subcommand(cmd_pgr::lav::make_subcommand())
        .subcommand(cmd_pgr::maf::make_subcommand())
        .subcommand(cmd_pgr::net::make_subcommand())
        .subcommand(cmd_pgr::psl::make_subcommand())
        .subcommand(cmd_pgr::pl::make_subcommand())
        .subcommand(cmd_pgr::twobit::make_subcommand())
        .subcommand(cmd_pgr::fa::make_subcommand())
        .subcommand(cmd_pgr::fas::make_subcommand())
        .subcommand(cmd_pgr::fq::make_subcommand())
        .after_help(
            r###"Subcommand groups:

* Sequences:
    * 2bit - Random access to .2bit files
    * fa   - FASTA operations: index, filter, stats
    * fas  - Block FASTA tools: consensus, variation
    * fq   - FASTQ operations: interleave, convert

* Genome alignments:
    * chain - Chain tools: sort, net, stitch
    * net   - Net tools: filter, syntenic
    * axt   - AXT conversion and sorting
    * lav   - LAV to PSL conversion
    * maf   - MAF to FASTA conversion
    * psl   - PSL operations: stats, to-chain

* Pipelines:
    * pl - Integrated pipelines: UCSC, TRF, etc.

"###,
        );

    // Check which subcomamnd the user ran...
    match app.get_matches().subcommand() {
        Some(("ms2dna", sub_matches)) => cmd_pgr::ms2dna::execute(sub_matches),
        Some(("axt", sub_matches)) => cmd_pgr::axt::execute(sub_matches),
        Some(("chain", sub_matches)) => cmd_pgr::chain::execute(sub_matches),
        Some(("lav", sub_matches)) => cmd_pgr::lav::execute(sub_matches),
        Some(("maf", sub_matches)) => cmd_pgr::maf::execute(sub_matches),
        Some(("net", sub_matches)) => cmd_pgr::net::execute(sub_matches),
        Some(("psl", sub_matches)) => cmd_pgr::psl::execute(sub_matches),
        Some(("pl", sub_matches)) => cmd_pgr::pl::execute(sub_matches),
        Some(("2bit", sub_matches)) => cmd_pgr::twobit::execute(sub_matches),
        Some(("fa", sub_matches)) => cmd_pgr::fa::execute(sub_matches),
        Some(("fas", sub_matches)) => cmd_pgr::fas::execute(sub_matches),
        Some(("fq", sub_matches)) => cmd_pgr::fq::execute(sub_matches),
        _ => unreachable!(),
    }?;

    Ok(())
}

// TODO: paralog
// TODO: fasr variation --indel
// TODO: fasr match
//  sparsemem -maxmatch -F -l %d -b -n -k 4 -threads 4 %s %s > %s
//  mummer -maxmatch -F -l %d -b -n %s %s > %s
//  $exe, $length, $genome, $query, $result
