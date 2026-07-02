//! Shared runtime helpers for `pgr paf` subcommands.
//!
//! Argument builders live in `crate::cmd_pgr::args`; this module provides
//! query execution and PAF output formatting shared across `paf` subcommands.

use clap::ArgMatches;

use pgr::libs::paf::fasta::{load_fasta_tsv, validate_tsv_covers_index, FastaStore};
use pgr::libs::paf::index::{PafIndex, QueryResult};
use pgr::libs::paf::msa_build::orient_interval;
use pgr::libs::paf::query::QueryOptions;

/// A region paired with its query results.
pub type QueryGroup = ((String, i32, i32), Vec<QueryResult>);

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
