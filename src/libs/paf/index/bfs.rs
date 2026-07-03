//! BFS transitive projection algorithm for [`super::PafIndex`].
//!
//! Walks the alignment graph outward from a seed region: at each depth, queries
//! both the forward index (`trees`) and the mirror index (`reverse_trees`) so
//! BFS can traverse from any sequence in both directions without requiring a
//! second PAF record with swapped roles.

use super::{project, PafIndex, PafMetadata, QueryResult};
use crate::libs::paf::cigar::gap_compressed_identity;
use coitrees::{Interval, IntervalNode, IntervalTree};
use std::collections::HashMap;

impl PafIndex {
    /// Transitive BFS: walks the alignment graph outward from a seed region up
    /// to `max_depth` hops, merging adjacent results within `merge_distance`.
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
            merge_results(&mut results, merge_distance);
        }
        results
    }
}

fn merge_results(results: &mut Vec<QueryResult>, max_gap: i32) {
    // Group by query_id, sort by query_start, merge adjacent within max_gap
    let mut groups: HashMap<u32, Vec<(usize, i32, i32)>> = HashMap::new();
    for (i, &(qid, q_iv, _t_iv, _, _, _, _)) in results.iter().enumerate() {
        groups
            .entry(qid)
            .or_default()
            .push((i, q_iv.first, q_iv.last));
    }
    let mut to_remove = Vec::new();
    for (_qid, mut items) in groups {
        if items.len() <= 1 {
            continue;
        }
        items.sort_by_key(|&(_, s, _)| s);
        let mut prev = items[0];
        for &curr in &items[1..] {
            // Merge when curr.start is within max_gap after prev.end (covers overlap + adjacency).
            if curr.1 <= prev.2 + max_gap {
                // Merge: keep the one with earlier target
                results[prev.0].1.last = results[prev.0].1.last.max(curr.2);
                prev.2 = prev.2.max(curr.2);
                to_remove.push(curr.0);
            } else {
                prev = curr;
            }
        }
    }
    to_remove.sort_unstable();
    to_remove.dedup();
    for &idx in to_remove.iter().rev() {
        results.remove(idx);
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
                55,
                55,
                '+',
            ),
        ];
        merge_results(&mut results, 10);
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
                100,
                100,
                '+',
            ),
        ];
        merge_results(&mut results, 10);
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
                50,
                50,
                '+',
            ),
        ];
        merge_results(&mut results, 10);
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
                55,
                55,
                '+',
            ),
            (
                0u32,
                Interval::new(65, 100, 0u32),
                Interval::new(65, 100, 1u32),
                vec![],
                65,
                65,
                '+',
            ),
        ];
        merge_results(&mut results, 10);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].1.first, 0);
        assert_eq!(results[0].1.last, 100);
    }
}
