use clap::*;
use pgr::libs::paf::index::QueryResult;

use super::query;

// Output BED3 (name start end), one line per query result.
fn output_bed(idx: &pgr::libs::paf::index::PafIndex, results: &[QueryResult]) {
    for (query_id, q_iv, _t_iv, _cigar, _, _) in results {
        let qname = idx.id_to_name(*query_id).unwrap_or("?");
        let (qs, qe) = if q_iv.first <= q_iv.last {
            (q_iv.first, q_iv.last)
        } else {
            (q_iv.last, q_iv.first)
        };
        println!("{qname}\t{qs}\t{qe}");
    }
}

pub fn make_subcommand() -> Command {
    query::add_query_args(Command::new("to-bed"))
        .about("Query PAF index and output BED3 coordinates")
        .after_help(
            r###"
Queries a PAF file or saved index (same logic as `pgr paf query`) and
outputs query coordinates as BED3 (name start end), one line per result.

This is the pipe-friendly coordinate-only view of `pgr paf query`.
All query options (region, --transitive, filters) are supported.

Notes:
* Input PAF files should contain cg:Z: tags for accurate projection
* Supports both plain text and gzipped (.gz) files (including BGZF)
* Reads from stdin if input file is 'stdin'

Examples:
1. Single region to BED:
   pgr paf to-bed alignments.paf chr1:1000-5000

2. Batch query from BED regions:
   pgr paf to-bed alignments.paf.idx -b regions.bed

3. With transitive BFS and identity filter:
   pgr paf to-bed alignments.paf chr1:1000-5000 -t --min-identity 0.8

"###,
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let (idx, all_results) = query::run_query(args)?;
    for (_, results) in &all_results {
        output_bed(&idx, results);
    }
    Ok(())
}
