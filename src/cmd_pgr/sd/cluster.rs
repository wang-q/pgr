//! `pgr sd cluster` — cluster overlapping SD hits and extract sequences.

use clap::{Arg, ArgMatches, Command};

/// Build the clap subcommand for cluster.
pub fn make_subcommand() -> Command {
    Command::new("cluster")
        .about("Clusters overlapping SD hits and extracts cluster FASTA")
        .after_help(
            r###"
Reads the refined SD PAF (from `pgr sd align`) and clusters overlapping
mates by union-find: both mates of one hit share a cluster, and intervals
overlapping on the same chromosome are unioned. Each cluster is written as a
FASTA file `cluster_N.fa` into the output directory, with headers in BISER
form `{species}#{chrom}{strand}#{start}#{end}` (0-based coordinates).

Examples:
1. Cluster refined SD hits:
   pgr sd cluster genome.fa hits.paf -o clusters.dir/
"###,
        )
        .arg(
            Arg::new("genome")
                .index(1)
                .required(true)
                .help("Genome FASTA file (plain or .gz)"),
        )
        .arg(
            Arg::new("paf")
                .index(2)
                .required(true)
                .help("Refined SD PAF file (from pgr sd align)"),
        )
        .arg(
            Arg::new("outdir")
                .long("outdir")
                .short('o')
                .default_value("clusters.dir")
                .help("Output directory for cluster FASTA files"),
        )
}
/// Execute the cluster command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let genome = args.get_one::<String>("genome").unwrap();
    let paf = args.get_one::<String>("paf").unwrap();
    let outdir = args.get_one::<String>("outdir").unwrap();

    let reader = pgr::reader(paf)?;
    let clusters = pgr::libs::sd::cluster::cluster_paf(reader, genome, outdir)?;
    log::info!("wrote {} cluster(s) to {}", clusters.len(), outdir);
    Ok(())
}
