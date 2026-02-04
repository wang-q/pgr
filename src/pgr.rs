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
        .subcommand(cmd_pgr::ms_to_dna::make_subcommand())
        .subcommand(cmd_pgr::axt::make_subcommand())
        .subcommand(cmd_pgr::chain::make_subcommand())
        .subcommand(cmd_pgr::chaining::make_subcommand())
        .subcommand(cmd_pgr::clust::make_subcommand())
        .subcommand(cmd_pgr::lav::make_subcommand())
        .subcommand(cmd_pgr::maf::make_subcommand())
        .subcommand(cmd_pgr::mat::make_subcommand())
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
    * 2bit - 2bit query and extraction
    * fa   - FASTA operations: info, records, transform, indexing
    * fas  - Block FA operations: info, subset, transform, file, variation
    * fq   - FASTQ interleaving and conversion

* Genome alignments:
    * chaining - Chaining alignments: psl
    * chain - Chain operations: sort, filter, transform, to-net
    * net   - Net operations: info, subset, transform, convert
    * axt   - AXT sorting and conversion
    * lav   - Convert to PSL
    * maf   - Convert to Block FA
    * psl   - PSL statistics, manipulation, and conversion

* Clustering:
    * clust - Algorithms: cc, dbscan, k-medoids, mcl

* Matrix:
    * mat   - Processing: compare, format, subset, to-pair, to-phylip

* Pipelines:
    * pl - Workflows: p2m, trf, ir, rept, ucsc

"###,
        );

    // Check which subcomamnd the user ran...
    match app.get_matches().subcommand() {
        Some(("ms-to-dna", sub_matches)) => cmd_pgr::ms_to_dna::execute(sub_matches),
        Some(("axt", sub_matches)) => cmd_pgr::axt::execute(sub_matches),
        Some(("chaining", sub_matches)) => cmd_pgr::chaining::execute(sub_matches),
        Some(("chain", sub_matches)) => cmd_pgr::chain::execute(sub_matches),
        Some(("clust", sub_matches)) => cmd_pgr::clust::execute(sub_matches),
        Some(("lav", sub_matches)) => cmd_pgr::lav::execute(sub_matches),
        Some(("maf", sub_matches)) => cmd_pgr::maf::execute(sub_matches),
        Some(("mat", sub_matches)) => cmd_pgr::mat::execute(sub_matches),
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
// TODO: fas variation --indel
// TODO: fas match
//  sparsemem -maxmatch -F -l %d -b -n -k 4 -threads 4 %s %s > %s
//  mummer -maxmatch -F -l %d -b -n %s %s > %s
//  $exe, $length, $genome, $query, $result
// TODO: Remove fully contained sequences
