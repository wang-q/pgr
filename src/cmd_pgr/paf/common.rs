use clap::{Arg, ArgMatches, Command};

use pgr::libs::paf::fasta::{load_fasta_tsv, validate_tsv_covers_index, FastaStore};
use pgr::libs::paf::index::{PafIndex, QueryResult};
use pgr::libs::poa::AlignmentParams;

use super::query;

/// A region paired with its query results.
pub type QueryGroup = ((String, i32, i32), Vec<QueryResult>);

/// Add common query arguments to a clap Command.
/// Shared by `paf query`, `paf to-bed`, `paf to-gfa`, `paf to-vcf`,
/// `paf to-fas`, and `paf to-maf`.
pub fn add_query_args(cmd: Command) -> Command {
    cmd.arg(
        Arg::new("infile")
            .required(true)
            .index(1)
            .help("Input PAF file or .paf.idx index to query"),
    )
    .arg(
        Arg::new("region")
            .index(2)
            .help("Target region to query (e.g. chr1:1000-5000)"),
    )
    .arg(
        Arg::new("bed_regions")
            .long("bed-regions")
            .short('b')
            .num_args(1)
            .help("BED file with multiple regions for batch query (name start end per line)"),
    )
    .arg(
        Arg::new("transitive")
            .long("transitive")
            .short('t')
            .num_args(0)
            .help("Enable transitive BFS traversal"),
    )
    .arg(
        Arg::new("max_depth")
            .long("max-depth")
            .short('m')
            .num_args(1)
            .default_value("2")
            .value_parser(clap::value_parser!(u16))
            .help("Maximum BFS depth (0 = unlimited, default: 2)"),
    )
    .arg(
        Arg::new("min_len")
            .long("min-len")
            .num_args(1)
            .default_value("10")
            .value_parser(clap::value_parser!(i32))
            .help("Minimum interval length to propagate (default: 10)"),
    )
    .arg(
        Arg::new("min_dist")
            .long("min-dist")
            .num_args(1)
            .default_value("10")
            .value_parser(clap::value_parser!(i32))
            .help("Minimum distance to merge adjacent intervals (default: 10)"),
    )
    .arg(
        Arg::new("min_identity")
            .long("min-identity")
            .num_args(1)
            .default_value("0.0")
            .value_parser(clap::value_parser!(f64))
            .help("Minimum gap-compressed identity (0.0-1.0, default: 0.0)"),
    )
    .arg(
        Arg::new("min_output_len")
            .long("min-output-len")
            .num_args(1)
            .default_value("0")
            .value_parser(clap::value_parser!(i32))
            .help("Minimum output interval length (default: 0 = no filter)"),
    )
    .arg(
        Arg::new("merge_distance")
            .long("merge-distance")
            .num_args(1)
            .default_value("0")
            .value_parser(clap::value_parser!(i32))
            .help("Merge adjacent output intervals within this distance (default: 0 = off)"),
    )
    .arg(
        Arg::new("min_degree")
            .long("min-degree")
            .num_args(1)
            .default_value("0")
            .value_parser(clap::value_parser!(usize))
            .help("Minimum distinct query sequences per region (default: 0 = off)"),
    )
    .arg(
        Arg::new("min_chain_length")
            .long("min-chain-length")
            .num_args(1)
            .default_value("0")
            .value_parser(clap::value_parser!(i32))
            .help("Minimum total aligned length per query (default: 0 = off)"),
    )
    .arg(
        Arg::new("subset_list")
            .long("subset-sequence-list")
            .num_args(1)
            .help("File with sequence names to include (one per line)"),
    )
    .arg(
        Arg::new("syntenic_filter")
            .long("syntenic-filter")
            .num_args(1)
            .help("UCSC chain file; drop query results whose query interval is not covered by any chain's query span (chain-level, both target and query name must match)"),
    )
}

/// Add POA scoring arguments to a clap Command.
/// Shared by `paf to-gfa`, `paf to-vcf`, `paf to-fas`, and `paf to-maf`.
pub fn add_poa_args(cmd: Command) -> Command {
    cmd.arg(
        Arg::new("match_score")
            .long("match")
            .num_args(1)
            .default_value("5")
            .value_parser(clap::value_parser!(i32))
            .allow_negative_numbers(true)
            .help("POA match score (default: 5)"),
    )
    .arg(
        Arg::new("mismatch_score")
            .long("mismatch")
            .num_args(1)
            .default_value("-4")
            .value_parser(clap::value_parser!(i32))
            .allow_negative_numbers(true)
            .help("POA mismatch score (default: -4)"),
    )
    .arg(
        Arg::new("gap_open")
            .long("gap-open")
            .num_args(1)
            .default_value("-8")
            .value_parser(clap::value_parser!(i32))
            .allow_negative_numbers(true)
            .help("POA gap open penalty (default: -8)"),
    )
    .arg(
        Arg::new("gap_extend")
            .long("gap-extend")
            .num_args(1)
            .default_value("-6")
            .value_parser(clap::value_parser!(i32))
            .allow_negative_numbers(true)
            .help("POA gap extend penalty (default: -6)"),
    )
}

/// Load fasta TSV, run query, validate TSV covers the index, build FastaStore.
/// Shared by `to-fas`, `to-gfa`, `to-maf`, `to-vcf`.
#[allow(clippy::type_complexity)]
pub fn prepare_query(args: &ArgMatches) -> anyhow::Result<(PafIndex, Vec<QueryGroup>, FastaStore)> {
    let tsv_path = args.get_one::<String>("fasta_tsv").unwrap();
    let seq_to_file = load_fasta_tsv(tsv_path)?;

    let (idx, all_results) = query::run_query(args)?;

    validate_tsv_covers_index(&seq_to_file, &idx)?;

    let fasta_store = FastaStore::new(&seq_to_file)?;

    Ok((idx, all_results, fasta_store))
}

/// Extract POA scoring parameters from ArgMatches.
pub fn get_poa_params(args: &ArgMatches) -> AlignmentParams {
    AlignmentParams {
        match_score: *args.get_one::<i32>("match_score").unwrap(),
        mismatch_score: *args.get_one::<i32>("mismatch_score").unwrap(),
        gap_open: *args.get_one::<i32>("gap_open").unwrap(),
        gap_extend: *args.get_one::<i32>("gap_extend").unwrap(),
    }
}
