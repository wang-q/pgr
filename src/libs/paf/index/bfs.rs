//! BFS transitive projection algorithm for [`super::PafIndex`].
//!
//! Walks the alignment graph outward from a seed region: at each depth, queries
//! both the forward index (`trees`) and the mirror index (`reverse_trees`) so
//! BFS can traverse from any sequence in both directions without requiring a
//! second PAF record with swapped roles.

use super::{project, PafIndex, PafMetadata, QueryResult};
use crate::libs::paf::cigar::{gap_compressed_identity, slice_cigar_by_target, CigarOp};
use crate::libs::paf::fasta::FastaStore;
use coitrees::{Interval, IntervalNode, IntervalTree};
use std::collections::HashMap;

impl PafIndex {
    /// Transitive BFS: walks the alignment graph outward from a seed region up
    /// to `max_depth` hops. Adjacent results are merged afterwards by the
    /// caller via [`PafIndex::merge_results`] when `--merge-distance` is set.
    #[allow(clippy::too_many_arguments)]
    pub fn query_transitive_bfs(
        &self,
        target_id: u32,
        start: i32,
        end: i32,
        max_depth: u16,
        min_len: usize,
        min_dist: i32,
        min_identity: f64,
        min_output_len: i32,
        merge_distance: i32,
        fasta_store: Option<&mut FastaStore>,
    ) -> Vec<QueryResult> {
        let mut results = Vec::new();
        let mut visited: HashMap<u32, SortedRanges> = HashMap::new();

        let init = visited
            .entry(target_id)
            .or_insert_with(|| SortedRanges::new(i32::MAX, min_dist))
            .insert(start, end);
        let mut current: Vec<(u32, i32, i32)> = init
            .into_iter()
            .filter(|&(s, e)| (e - s).unsigned_abs() >= min_len as u32)
            .map(|(s, e)| (target_id, s, e))
            .collect();

        let mut depth = 0u16;
        while !current.is_empty() && (max_depth == 0 || depth < max_depth) {
            let mut next = vec![];
            for &(tid, cs, ce) in &current {
                // Query both the forward index (`trees`) and the mirror index
                // (`reverse_trees`) so BFS can traverse from any sequence in
                // both directions. Arc clones are cheap (refcount bump) and
                // release the borrow on `self` before the closure captures it.
                let fwd = self.trees.get(&tid).cloned();
                let rev = self.reverse_trees.get(&tid).cloned();
                for tree in fwd.into_iter().chain(rev) {
                    tree.query(cs, ce, |iv: &IntervalNode<PafMetadata, u32>| {
                        let m = &iv.metadata;
                        let os = cs.max(iv.first);
                        let oe = ce.min(iv.last);
                        if os >= oe {
                            return;
                        }
                        let cigar = self.resolve_cigar(&m.cigar);
                        if gap_compressed_identity(&cigar) < min_identity {
                            return;
                        }
                        if let Some((qs, qe, ts, te)) = project(os, oe, m, &cigar) {
                            if (qe - qs).abs() < min_output_len {
                                return;
                            }
                            results.push((
                                m.query_id,
                                Interval::new(qs, qe, m.query_id),
                                Interval::new(ts, te, tid),
                                cigar,
                                m.target_start,
                                m.query_start,
                                m.strand,
                            ));
                            if m.query_id != tid {
                                let sr = visited
                                    .entry(m.query_id)
                                    .or_insert_with(|| SortedRanges::new(i32::MAX, min_dist));
                                for (ns, ne) in sr.insert(qs, qe) {
                                    if (ne - ns).unsigned_abs() >= min_len as u32 {
                                        next.push((m.query_id, ns, ne));
                                    }
                                }
                            }
                        }
                    });
                }
            }
            current = next;
            depth += 1;
        }
        if merge_distance > 0 {
            self.merge_results(&mut results, merge_distance, fasta_store);
        }
        results
    }
}

impl PafIndex {
    fn merge_results(
        &self,
        results: &mut Vec<QueryResult>,
        max_gap: i32,
        mut fasta_store: Option<&mut FastaStore>,
    ) {
        if results.len() < 2 {
            return;
        }
        // Group by (query_id, target_id, strand). Within a group, adjacent
        // results are merged; the original-record check uses rec_ts/rec_qs.
        let mut groups: HashMap<(u32, u32, char), Vec<usize>> = HashMap::new();
        for (i, (qid, _q_iv, t_iv, _cigar, _rec_ts, _rec_qs, strand)) in results.iter().enumerate()
        {
            groups
                .entry((*qid, t_iv.metadata, *strand))
                .or_default()
                .push(i);
        }

        let mut to_remove = Vec::new();
        for (_key, mut idxs) in groups {
            if idxs.len() <= 1 {
                continue;
            }
            idxs.sort_by_key(|&i| results[i].1.first);
            let mut curr = idxs[0];
            for &next in &idxs[1..] {
                let (c_qs, c_qe) = (results[curr].1.first, results[curr].1.last);
                let (n_qs, n_qe) = (results[next].1.first, results[next].1.last);
                let (c_ts, c_te) = (results[curr].2.first, results[curr].2.last);
                let (n_ts, n_te) = (results[next].2.first, results[next].2.last);

                if n_qs > c_qe + max_gap {
                    curr = next;
                    continue;
                }

                let merged_qs = c_qs.min(n_qs);
                let merged_qe = c_qe.max(n_qe);
                let merged_ts = c_ts.min(n_ts);
                let merged_te = c_te.max(n_te);

                let same_record = results[curr].4 == results[next].4
                    && results[curr].5 == results[next].5
                    && results[curr].6 == results[next].6;

                let new_cigar = if same_record {
                    Some(slice_cigar_by_target(
                        &results[curr].3,
                        results[curr].4,
                        merged_ts,
                        merged_te,
                    ))
                } else {
                    fasta_store.as_mut().and_then(|store| {
                        Self::recompute_gap_cigar(
                            store,
                            self.id_to_name(results[curr].0).unwrap_or("?"),
                            self.id_to_name(results[curr].2.metadata).unwrap_or("?"),
                            merged_qs,
                            merged_qe,
                            merged_ts,
                            merged_te,
                            results[curr].6,
                        )
                    })
                };

                if let Some(cigar) = new_cigar {
                    results[curr].1.first = merged_qs;
                    results[curr].1.last = merged_qe;
                    results[curr].2.first = merged_ts;
                    results[curr].2.last = merged_te;
                    results[curr].3 = cigar;
                    to_remove.push(next);
                } else {
                    curr = next;
                }
            }
        }

        to_remove.sort_unstable();
        to_remove.dedup();
        for &idx in to_remove.iter().rev() {
            results.remove(idx);
        }
    }

    /// Recompute a gap-filled CIGAR between two merged result intervals.
    ///
    /// Fetches the merged target/query ranges from `store` and emits `=`/`X`
    /// ops. Returns `None` if the sequences cannot be fetched or if the
    /// merged ranges contain indels (unequal lengths), in which case merging
    /// is conservatively skipped.
    #[allow(clippy::too_many_arguments)]
    fn recompute_gap_cigar(
        store: &mut FastaStore,
        qname: &str,
        tname: &str,
        qs: i32,
        qe: i32,
        ts: i32,
        te: i32,
        strand: char,
    ) -> Option<Vec<CigarOp>> {
        let (t_seq, _) = store.fetch_range(tname, ts, te).ok()?;
        let (q_seq_fwd, _) = store.fetch_range(qname, qs, qe).ok()?;
        let q_seq: Vec<u8> = if strand == '-' {
            crate::libs::nt::rev_comp(&q_seq_fwd).collect()
        } else {
            q_seq_fwd
        };
        if t_seq.len() != q_seq.len() {
            return None;
        }
        let mut ops: Vec<CigarOp> = Vec::new();
        for (&t, &q) in t_seq.iter().zip(q_seq.iter()) {
            let op = if t.eq_ignore_ascii_case(&q) { '=' } else { 'X' };
            match ops.last_mut() {
                Some(last) if last.op() == op => {
                    *last = CigarOp::new(last.len() + 1, op);
                }
                _ => ops.push(CigarOp::new(1, op)),
            }
        }
        Some(ops)
    }
}

struct SortedRanges {
    ranges: Vec<(i32, i32)>,
    seq_len: i32,
    min_dist: i32,
}

impl SortedRanges {
    fn new(seq_len: i32, min_dist: i32) -> Self {
        Self {
            ranges: vec![],
            seq_len,
            min_dist,
        }
    }

    fn insert(&mut self, start: i32, end: i32) -> Vec<(i32, i32)> {
        let (mut s, mut e) = if start <= end {
            (start, end)
        } else {
            (end, start)
        };
        s = s.max(0);
        e = e.min(self.seq_len);
        if s >= e {
            return vec![];
        }
        // Find insertion point
        let mut i = match self.ranges.binary_search_by_key(&s, |&(a, _)| a) {
            Ok(p) => p,
            Err(p) => p,
        };
        // Step back if previous range overlaps
        if i > 0 && self.ranges[i - 1].1 > s {
            i -= 1;
        }
        // Expand to nearby ranges
        if i > 0 && (s - self.ranges[i - 1].1).abs() <= self.min_dist {
            s = self.ranges[i - 1].1;
            i -= 1;
        }
        if i < self.ranges.len() && (self.ranges[i].0 - e).abs() <= self.min_dist {
            e = self.ranges[i].0;
        }
        // Collect non-covered portions
        let mut out = vec![];
        let mut cur = s;
        let mut j = i;
        while j < self.ranges.len() && cur < e {
            let (rs, re) = self.ranges[j];
            if rs > e {
                break;
            }
            if cur < rs {
                out.push((cur, rs.min(e)));
            }
            cur = cur.max(re);
            j += 1;
        }
        if cur < e {
            out.push((cur, e));
        }
        // Merge
        let first = match self.ranges.binary_search_by_key(&s, &|&(a, _)| a) {
            Ok(p) => p,
            Err(p) => {
                if p > 0 && self.ranges[p - 1].1 >= s {
                    p - 1
                } else {
                    p
                }
            }
        };
        let last = self.ranges[first..]
            .iter()
            .position(|&(a, _)| a > e)
            .map(|p| first + p)
            .unwrap_or(self.ranges.len());
        let old_end = self.ranges.get(last.saturating_sub(1)).map(|r| r.1);
        let merged = (
            s.min(self.ranges.get(first).map(|r| r.0).unwrap_or(s)),
            e.max(old_end.unwrap_or(e)),
        );
        self.ranges.drain(first..last);
        self.ranges.insert(first, merged);
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use coitrees::Interval;
    use indexmap::IndexMap;

    fn empty_index() -> PafIndex {
        PafIndex {
            names: IndexMap::new(),
            trees: HashMap::new(),
            reverse_trees: HashMap::new(),
            lazy_source: None,
            lazy_source_path: None,
        }
    }

    fn tmp_fasta_store(seqs: &[(&str, &str)]) -> (tempfile::TempDir, FastaStore) {
        use std::io::Write;
        let dir = tempfile::tempdir().unwrap();
        let mut entries = IndexMap::new();
        for (name, seq) in seqs {
            let path = dir.path().join(format!("{name}.fa.gz"));
            let file = std::fs::File::create(&path).unwrap();
            let mut writer = noodles_bgzf::io::Writer::new(file);
            writeln!(writer, ">{name}").unwrap();
            writeln!(writer, "{seq}").unwrap();
            writer.flush().unwrap();
            drop(writer);
            crate::libs::fmt::fa::build_gzi_index(path.to_str().unwrap()).unwrap();
            entries.insert(name.to_string(), path.to_string_lossy().to_string());
        }
        let store = FastaStore::new(&entries).unwrap();
        (dir, store)
    }

    #[test]
    fn test_sorted_ranges_disjoint() {
        let mut sr = SortedRanges::new(1000, 0);
        assert_eq!(sr.insert(100, 200), vec![(100, 200)]);
        assert_eq!(sr.insert(300, 400), vec![(300, 400)]);
    }

    #[test]
    fn test_sorted_ranges_contained() {
        let mut sr = SortedRanges::new(1000, 0);
        sr.insert(100, 200);
        assert!(sr.insert(120, 180).is_empty());
    }

    #[test]
    fn test_sorted_ranges_adjacent() {
        let mut sr = SortedRanges::new(1000, 10);
        sr.insert(100, 200);
        sr.insert(210, 300);
        assert_eq!(sr.ranges.len(), 1);
        assert_eq!(sr.ranges[0], (100, 300));
    }

    #[test]
    fn test_merge_adjacent_intervals() {
        let mut results = vec![
            (
                0u32,
                Interval::new(0, 50, 0u32),
                Interval::new(0, 50, 1u32),
                vec![],
                0,
                0,
                '+',
            ),
            (
                0u32,
                Interval::new(55, 100, 0u32),
                Interval::new(55, 100, 1u32),
                vec![],
                0,
                0,
                '+',
            ),
        ];
        empty_index().merge_results(&mut results, 10, None);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].1.first, 0);
        assert_eq!(results[0].1.last, 100);
    }

    #[test]
    fn test_merge_no_merge_when_far() {
        let mut results = vec![
            (
                0u32,
                Interval::new(0, 50, 0u32),
                Interval::new(0, 50, 1u32),
                vec![],
                0,
                0,
                '+',
            ),
            (
                0u32,
                Interval::new(100, 150, 0u32),
                Interval::new(100, 150, 1u32),
                vec![],
                0,
                0,
                '+',
            ),
        ];
        empty_index().merge_results(&mut results, 10, None);
        assert_eq!(results.len(), 2);
    }

    #[test]
    fn test_merge_overlapping_intervals() {
        // Overlapping intervals must merge regardless of max_gap (regression test
        // for abs() bug that wrongly rejected overlaps when overlap > max_gap).
        let mut results = vec![
            (
                0u32,
                Interval::new(0, 100, 0u32),
                Interval::new(0, 100, 1u32),
                vec![],
                0,
                0,
                '+',
            ),
            (
                0u32,
                Interval::new(50, 150, 0u32),
                Interval::new(50, 150, 1u32),
                vec![],
                0,
                0,
                '+',
            ),
        ];
        empty_index().merge_results(&mut results, 10, None);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].1.first, 0);
        assert_eq!(results[0].1.last, 150);
    }

    #[test]
    fn test_merge_chained_intervals() {
        // Three intervals where each pair is within max_gap, but the third is
        // only mergeable after the second extends prev.end (regression test for
        // prev.2 not being updated after merge).
        let mut results = vec![
            (
                0u32,
                Interval::new(0, 50, 0u32),
                Interval::new(0, 50, 1u32),
                vec![],
                0,
                0,
                '+',
            ),
            (
                0u32,
                Interval::new(55, 60, 0u32),
                Interval::new(55, 60, 1u32),
                vec![],
                0,
                0,
                '+',
            ),
            (
                0u32,
                Interval::new(65, 100, 0u32),
                Interval::new(65, 100, 1u32),
                vec![],
                0,
                0,
                '+',
            ),
        ];
        empty_index().merge_results(&mut results, 10, None);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].1.first, 0);
        assert_eq!(results[0].1.last, 100);
    }

    #[test]
    fn test_merge_same_record_cigar_recomputed() {
        // Same source record: merging must rebuild the CIGAR for the union of
        // target intervals, not just extend query coordinates.
        let cigar = crate::libs::paf::cigar::parse_cigar("25=5X25=").unwrap();
        let mut results = vec![
            (
                0u32,
                Interval::new(0, 25, 0u32),
                Interval::new(0, 25, 1u32),
                cigar.clone(),
                0,
                0,
                '+',
            ),
            (
                0u32,
                Interval::new(30, 55, 0u32),
                Interval::new(30, 55, 1u32),
                cigar.clone(),
                0,
                0,
                '+',
            ),
        ];
        empty_index().merge_results(&mut results, 10, None);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].1.first, 0);
        assert_eq!(results[0].1.last, 55);
        assert_eq!(results[0].2.first, 0);
        assert_eq!(results[0].2.last, 55);
        let out = crate::libs::paf::cigar::format_cigar(&results[0].3);
        assert_eq!(out, "25=5X25=");
    }

    #[test]
    fn test_merge_without_fasta_does_not_merge_different_records() {
        // Different source records without a FastaStore: merging must be
        // skipped because the gap CIGAR cannot be recomputed.
        let mut results = vec![
            (
                0u32,
                Interval::new(0, 25, 0u32),
                Interval::new(0, 25, 1u32),
                vec![],
                0,
                0,
                '+',
            ),
            (
                0u32,
                Interval::new(30, 55, 0u32),
                Interval::new(30, 55, 1u32),
                vec![],
                30,
                0,
                '+',
            ),
        ];
        empty_index().merge_results(&mut results, 10, None);
        assert_eq!(results.len(), 2);
    }

    #[test]
    fn test_merge_different_records_with_fasta() {
        // Different source records with a FastaStore: the gap is filled by
        // comparing the merged target/query sequences base-by-base.
        let (_dir, mut store) = tmp_fasta_store(&[
            ("T", "AAAAAAAAAACCCCCCCCCCTTTTTTTTTT"),
            ("Q", "AAAAAAAAAACCCCCCCCCgTTTTTTTTTT"),
        ]);
        let mut idx = empty_index();
        idx.names.insert("Q".to_string(), 0);
        idx.names.insert("T".to_string(), 1);

        let mut results = vec![
            (
                0u32,
                Interval::new(0, 20, 0u32),
                Interval::new(0, 20, 1u32),
                vec![],
                0,
                0,
                '+',
            ),
            (
                0u32,
                Interval::new(20, 30, 0u32),
                Interval::new(20, 30, 1u32),
                vec![],
                20,
                20,
                '+',
            ),
        ];
        idx.merge_results(&mut results, 10, Some(&mut store));
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].1.first, 0);
        assert_eq!(results[0].1.last, 30);
        assert_eq!(results[0].2.first, 0);
        assert_eq!(results[0].2.last, 30);
        let out = crate::libs::paf::cigar::format_cigar(&results[0].3);
        assert_eq!(out, "19=1X10=");
    }

    #[test]
    fn test_merge_different_records_with_fasta_minus_strand() {
        // Different source records on the '-' strand: query coords are forward,
        // but the merged CIGAR must compare the reverse complement of Q to T.
        let (_dir, mut store) = tmp_fasta_store(&[("T", "AAAACCCCGGGG"), ("Q", "CCCCGGGGTTTT")]);
        let mut idx = empty_index();
        idx.names.insert("Q".to_string(), 0);
        idx.names.insert("T".to_string(), 1);

        let mut results = vec![
            (
                0u32,
                Interval::new(0, 6, 0u32),
                Interval::new(0, 6, 1u32),
                vec![],
                0,
                0,
                '-',
            ),
            (
                0u32,
                Interval::new(6, 12, 0u32),
                Interval::new(6, 12, 1u32),
                vec![],
                6,
                6,
                '-',
            ),
        ];
        idx.merge_results(&mut results, 10, Some(&mut store));
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].1.first, 0);
        assert_eq!(results[0].1.last, 12);
        assert_eq!(results[0].2.first, 0);
        assert_eq!(results[0].2.last, 12);
        let out = crate::libs::paf::cigar::format_cigar(&results[0].3);
        assert_eq!(out, "12=");
    }
}
