//! Query engine and CIGAR resolution for [`super::PafIndex`].

use super::{CigarStore, PafIndex, PafMetadata};
use crate::libs::paf::cigar::{extract_cigar, gap_compressed_identity, reverse_cigar, CigarOp};
use crate::libs::paf::parser::parse_paf_line;
use coitrees::{Interval, IntervalNode, IntervalTree};
use noodles_bgzf as bgzf;
use std::io::BufRead;

impl PafIndex {
    /// Fetch CIGAR ops for a lazy record by seeking to its BGZF virtual position.
    pub fn fetch_cigar(&self, vpos: u64) -> Vec<CigarOp> {
        if let Some(ref src) = self.lazy_source {
            // Recover from a poisoned mutex rather than panicking (Zero Panic).
            let mut reader = src.lock().unwrap_or_else(|e| e.into_inner());
            if reader.seek(bgzf::VirtualPosition::from(vpos)).is_err() {
                return vec![];
            }
            let mut line = String::new();
            if reader.read_line(&mut line).is_err() {
                return vec![];
            }
            // Parse tags from the line and extract CIGAR.
            let trimmed = line.trim_end_matches('\n');
            if let Ok(rec) = parse_paf_line(trimmed) {
                return extract_cigar(&rec.tags).unwrap_or_default();
            }
            vec![]
        } else {
            vec![]
        }
    }

    /// Resolve a `CigarStore` to owned CIGAR ops (clone if owned, fetch if lazy).
    pub(crate) fn resolve_cigar(&self, store: &CigarStore) -> Vec<CigarOp> {
        match store {
            CigarStore::Owned(c) => c.clone(),
            CigarStore::Lazy(vpos) => self.fetch_cigar(*vpos),
            CigarStore::LazyReversed(vpos) => reverse_cigar(&self.fetch_cigar(*vpos)),
        }
    }

    /// Single-hop query: returns all records overlapping `[start, end)` on
    /// `target_id`, filtered by `min_identity` and `min_output_len`.
    pub fn query(
        &self,
        target_id: u32,
        start: i32,
        end: i32,
        min_identity: f64,
        min_output_len: i32,
    ) -> Vec<super::QueryResult> {
        let mut results = Vec::new();
        if let Some(tree) = self.trees.get(&target_id) {
            tree.query(start, end, |iv: &IntervalNode<PafMetadata, u32>| {
                let m = &iv.metadata;
                let cigar = self.resolve_cigar(&m.cigar);
                if gap_compressed_identity(&cigar) < min_identity {
                    return;
                }
                if let Some((qs, qe, ts, te)) = project(start, end, m, &cigar) {
                    if (qe - qs).abs() < min_output_len {
                        return;
                    }
                    results.push((
                        m.query_id,
                        Interval::new(qs, qe, m.query_id),
                        Interval::new(ts, te, target_id),
                        cigar,
                        m.target_start,
                        m.query_start,
                        m.strand,
                    ));
                }
            });
        }
        results
    }
}

/// Project a target sub-interval `[ts, te)` onto the query coordinates.
///
/// Returns `(q_start, q_end, t_start, t_end)` of the union of overlapping
/// CIGAR op spans, or `None` if no overlap.
pub(crate) fn project(
    ts: i32,
    te: i32,
    m: &PafMetadata,
    cigar: &[CigarOp],
) -> Option<(i32, i32, i32, i32)> {
    if cigar.is_empty() {
        let off = (ts - m.target_start).max(0);
        let len = (te - ts).min((m.target_end - m.target_start) - off);
        if len <= 0 {
            return None;
        }
        return Some((
            m.query_start + off,
            m.query_start + off + len,
            ts.max(m.target_start),
            (ts + len).min(m.target_end),
        ));
    }
    // Walk all CIGAR ops; accumulate the union of query/target intervals that
    // overlap [ts, te). Returning at the first overlap would truncate the
    // projection when the query spans indels (e.g. 4=3I3= over [0,7) must
    // project to query [0,10), not [0,4)).
    //
    // `cq` is the CIGAR query offset (0-based from the start of the aligned
    // query region). For '+' strand records this maps directly to forward
    // query coordinates (`query_start + cq`). For '-' strand records the
    // CIGAR describes RC(query) vs target, so RC offset `[rc_lo, rc_hi)` maps
    // to forward query `[query_end - rc_hi, query_end - rc_lo)`. PAF stores
    // `query_start`/`query_end` as forward-strand coordinates regardless of
    // strand, so this conversion is needed for sub-interval queries.
    let mut ct = m.target_start;
    let mut cq: i32 = 0;
    let mut q_min = i32::MAX;
    let mut q_max = i32::MIN;
    let mut t_min = i32::MAX;
    let mut t_max = i32::MIN;
    let mut found = false;
    // Convert RC offset interval [rc_lo, rc_hi) to forward query coordinates.
    let rc_to_forward = |rc_lo: i32, rc_hi: i32| -> (i32, i32) {
        if m.strand == '-' {
            (m.query_end - rc_hi, m.query_end - rc_lo)
        } else {
            (m.query_start + rc_lo, m.query_start + rc_hi)
        }
    };
    for op in cigar {
        let td = op.target_delta() as i32;
        let qd = op.query_delta() as i32;
        let ss = ct;
        let se = ct + td;
        match op.op() {
            '=' | 'X' | 'M' => {
                let os = ts.max(ss);
                let oe = te.min(se);
                if os < oe {
                    let off = os - ss;
                    let len = oe - os;
                    let (qs, qe) = rc_to_forward(cq + off, cq + off + len);
                    q_min = q_min.min(qs);
                    q_max = q_max.max(qe);
                    t_min = t_min.min(os);
                    t_max = t_max.max(oe);
                    found = true;
                }
            }
            'I' => {
                // Insertion in query at target position ct.
                // Include when ct lies within the queried target span.
                if ct >= ts && ct < te {
                    let (qs, qe) = rc_to_forward(cq, cq + qd);
                    q_min = q_min.min(qs);
                    q_max = q_max.max(qe);
                    found = true;
                }
            }
            'D' => {
                let os = ts.max(ss);
                let oe = te.min(se);
                if os < oe {
                    t_min = t_min.min(os);
                    t_max = t_max.max(oe);
                    found = true;
                }
            }
            _ => {}
        }
        ct = se;
        cq += qd;
    }
    if found {
        Some((q_min, q_max, t_min, t_max))
    } else {
        None
    }
}
