//! BED3 output helpers for PAF query results.

use std::io::Write;

use super::index::{PafIndex, QueryResult};
use super::msa_build::orient_interval;

/// Write BED3 (`name start end`), one line per query result.
///
/// Missing query ids are logged and emitted as `"?"` to keep output flowing.
pub fn write_bed3(
    idx: &PafIndex,
    results: &[QueryResult],
    writer: &mut dyn Write,
) -> anyhow::Result<()> {
    for (query_id, q_iv, _t_iv, _cigar, _, _, _) in results {
        let qname = idx.id_to_name(*query_id).unwrap_or_else(|| {
            log::warn!("query id {} not found in index", query_id);
            "?"
        });
        let (qs, qe) = orient_interval(q_iv.first, q_iv.last);
        writeln!(writer, "{qname}\t{qs}\t{qe}")?;
    }
    Ok(())
}
