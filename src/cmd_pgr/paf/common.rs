use clap::ArgMatches;

use pgr::libs::paf::fasta::{load_fasta_tsv, validate_tsv_covers_index, FastaStore};
use pgr::libs::paf::index::{PafIndex, QueryResult};
use pgr::libs::poa::AlignmentParams;

use super::query;

/// A region paired with its query results.
pub type QueryGroup = ((String, i32, i32), Vec<QueryResult>);

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
