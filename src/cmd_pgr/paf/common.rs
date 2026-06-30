use clap::{Arg, ArgMatches, Command};

use pgr::libs::paf::fasta::{load_fasta_tsv, validate_tsv_covers_index, FastaStore};
use pgr::libs::paf::index::{PafIndex, QueryResult};
use pgr::libs::paf::msa_build::orient_interval;
use pgr::libs::paf::query::QueryOptions;

// Re-export shared POA argument builder and extractor from `cmd_pgr::args`.
pub use crate::cmd_pgr::args::{add_poa_args, get_poa_params};

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

/// Add the required `-f/--fasta-tsv` argument.
/// Shared by `paf to-gfa`, `paf to-vcf`, `paf to-fas`, and `paf to-maf`.
pub fn add_fasta_tsv_arg(cmd: Command) -> Command {
    cmd.arg(
        Arg::new("fasta_tsv")
            .long("fasta-tsv")
            .short('f')
            .required(true)
            .num_args(1)
            .help("TSV file: genome_name <tab> bgzf_fasta_path"),
    )
}

/// Add the `--msa` flag for POA-based multi-way output.
/// Shared by `paf to-fas` and `paf to-maf`.
pub fn add_msa_flag(cmd: Command) -> Command {
    cmd.arg(
        Arg::new("msa")
            .long("msa")
            .num_args(0)
            .help("Merge results per region into a multi-way block via POA"),
    )
}

/// Extract [`QueryOptions`] from clap matches.
pub fn query_options_from_args(args: &ArgMatches) -> QueryOptions {
    QueryOptions {
        infile: args.get_one::<String>("infile").unwrap().clone(),
        region: args.get_one::<String>("region").cloned(),
        bed_regions: args.get_one::<String>("bed_regions").cloned(),
        transitive: args.get_flag("transitive"),
        max_depth: *args.get_one::<u16>("max_depth").unwrap(),
        min_len: *args.get_one::<i32>("min_len").unwrap(),
        min_dist: *args.get_one::<i32>("min_dist").unwrap(),
        min_identity: *args.get_one::<f64>("min_identity").unwrap(),
        min_output_len: *args.get_one::<i32>("min_output_len").unwrap(),
        merge_distance: *args.get_one::<i32>("merge_distance").unwrap(),
        min_degree: *args.get_one::<usize>("min_degree").unwrap(),
        min_chain_length: *args.get_one::<i32>("min_chain_length").unwrap(),
        subset_list: args.get_one::<String>("subset_list").cloned(),
        syntenic_filter: args.get_one::<String>("syntenic_filter").cloned(),
    }
}

/// Shared query logic: parse args, build/load index, run queries, apply filters.
/// Returns the index and a list of (region, results) pairs.
#[allow(clippy::type_complexity)]
pub fn run_query(
    args: &ArgMatches,
) -> anyhow::Result<(PafIndex, Vec<((String, i32, i32), Vec<QueryResult>)>)> {
    let opts = query_options_from_args(args);
    pgr::libs::paf::query::run_query(&opts)
}

// Output PAF: 12 columns + gi/bi/cg tags.
pub fn output_paf(idx: &PafIndex, results: &[QueryResult]) {
    for (query_id, q_iv, t_iv, cigar, _, _, strand) in results {
        let qname = idx.id_to_name(*query_id).unwrap_or("?");
        let tname = idx.id_to_name(t_iv.metadata).unwrap_or("?");
        let block_len = (q_iv.last - q_iv.first).abs().max(1) as u32;
        let matches = pgr::libs::paf::cigar::cigar_stats(cigar).matches;
        let gi = pgr::libs::paf::cigar::gap_compressed_identity(cigar);
        let bi = pgr::libs::paf::cigar::block_identity(cigar);
        let cg = pgr::libs::paf::cigar::format_cigar(cigar);
        let (qs, qe) = orient_interval(q_iv.first, q_iv.last);
        let (ts, te) = orient_interval(t_iv.first, t_iv.last);
        println!(
            "{}\t0\t{}\t{}\t{}\t{}\t0\t{}\t{}\t{}\t{}\t255\tgi:f:{:.6}\tbi:f:{:.6}\tcg:Z:{}",
            qname, qs, qe, strand, tname, ts, te, matches, block_len, gi, bi, cg
        );
    }
}

/// Load fasta TSV, run query, validate TSV covers the index, build FastaStore.
/// Shared by `to-fas`, `to-gfa`, `to-maf`, `to-vcf`.
#[allow(clippy::type_complexity)]
pub fn prepare_query(args: &ArgMatches) -> anyhow::Result<(PafIndex, Vec<QueryGroup>, FastaStore)> {
    let tsv_path = args.get_one::<String>("fasta_tsv").unwrap();
    let seq_to_file = load_fasta_tsv(tsv_path)?;

    let (idx, all_results) = run_query(args)?;

    validate_tsv_covers_index(&seq_to_file, &idx)?;

    let fasta_store = FastaStore::new(&seq_to_file)?;

    Ok((idx, all_results, fasta_store))
}
