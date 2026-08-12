//! Greedy best-edge layout of the overlap graph (OLC stage S2).
//!
//! Each unitig has two ends (0 = position 0, 1 = position len). Overlap
//! records are normalized into directed extension edges between ends with a
//! flip bit (target read reverse relative to the source). Layouts are built
//! by seeding unplaced unitigs longest-first and extending both directions
//! through mutual-best edges, stopping at placed unitigs, ambiguous
//! junctions (repeat ends), and non-reciprocal edges.

use super::overlap::{Overlap, Unitig};

/// One unitig on a layout path.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LayoutStep {
    /// Unitig index.
    pub unitig: usize,
    /// Orientation of the unitig in the contig (`+` / `-`).
    pub strand: char,
    /// Interval of the unitig in the contig (0-based, half-open).
    pub q_start: usize,
    pub q_end: usize,
    /// Exact overlap with the previous step (0 for the first step).
    pub overlap_len: usize,
}

/// An ordered layout path of unitigs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Layout {
    pub steps: Vec<LayoutStep>,
}

/// A directed extension edge between two unitig ends.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Edge {
    from: usize,
    from_end: usize,
    to: usize,
    to_end: usize,
    flip: bool,
    length: usize,
}

/// Builds greedy layouts from verified overlaps.
///
/// Seeds are processed longest-first (length then id, deterministic).
/// Extension requires a mutual-best junction: the target's best edge at the
/// junction end must point back to the current unitig. A unitig end whose
/// two best edges are near-equal and go to different targets is treated as
/// a repeat junction and blocks extension. Output layouts are sorted by
/// total length descending.
pub fn build_layouts(unitigs: &[Unitig], overlaps: &[Overlap]) -> anyhow::Result<Vec<Layout>> {
    let n = unitigs.len();
    let mut ends: Vec<Vec<Edge>> = vec![Vec::new(); 2 * n];
    for ov in overlaps {
        let qlen = unitigs[ov.qid].seq.len();
        let tlen = unitigs[ov.tid].seq.len();
        let flip = ov.strand == '-';
        if flip {
            // q.B <-> t.B (q suffix matches rc(t) prefix).
            if ov.q_end == qlen && ov.t_end == tlen {
                push_edge(&mut ends, ov.qid, 1, ov.tid, 1, true, ov.length);
                push_edge(&mut ends, ov.tid, 1, ov.qid, 1, true, ov.length);
            }
            // q.A <-> t.A (q prefix matches rc(t) suffix).
            if ov.q_start == 0 && ov.t_start == 0 {
                push_edge(&mut ends, ov.qid, 0, ov.tid, 0, true, ov.length);
                push_edge(&mut ends, ov.tid, 0, ov.qid, 0, true, ov.length);
            }
        } else {
            // q.B <-> t.A (q suffix matches t prefix).
            if ov.q_end == qlen && ov.t_start == 0 {
                push_edge(&mut ends, ov.qid, 1, ov.tid, 0, false, ov.length);
                push_edge(&mut ends, ov.tid, 0, ov.qid, 1, false, ov.length);
            }
            // q.A <-> t.B (q prefix matches t suffix).
            if ov.q_start == 0 && ov.t_end == tlen {
                push_edge(&mut ends, ov.qid, 0, ov.tid, 1, false, ov.length);
                push_edge(&mut ends, ov.tid, 1, ov.qid, 0, false, ov.length);
            }
        }
    }
    for e in &mut ends {
        dedup_and_sort(e);
    }

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        unitigs[b]
            .seq
            .len()
            .cmp(&unitigs[a].seq.len())
            .then(a.cmp(&b))
    });
    let mut placed = vec![false; n];
    let mut layouts = Vec::new();
    for &seed in &order {
        if placed[seed] {
            continue;
        }
        placed[seed] = true;
        let mut steps = vec![LayoutStep {
            unitig: seed,
            strand: '+',
            q_start: 0,
            q_end: unitigs[seed].seq.len(),
            overlap_len: 0,
        }];
        extend(&mut steps, &ends, &mut placed, seed, 1, '+', true);
        extend(&mut steps, &ends, &mut placed, seed, 0, '+', false);
        layouts.push(Layout { steps });
    }

    for layout in &mut layouts {
        let mut pos = 0usize;
        for (i, step) in layout.steps.iter_mut().enumerate() {
            if i == 0 {
                step.q_start = 0;
                step.q_end = unitigs[step.unitig].seq.len();
                pos = step.q_end;
            } else {
                let prev_end = pos;
                anyhow::ensure!(
                    step.overlap_len <= prev_end,
                    "layout {}: overlap {} exceeds the previous contig end {}",
                    i,
                    step.overlap_len,
                    prev_end
                );
                step.q_start = prev_end - step.overlap_len;
                step.q_end = step.q_start + unitigs[step.unitig].seq.len();
                pos = step.q_end;
            }
        }
    }
    layouts.sort_by(|a, b| {
        let la: usize = a.steps.iter().map(|s| s.q_end - s.q_start).sum();
        let lb: usize = b.steps.iter().map(|s| s.q_end - s.q_start).sum();
        lb.cmp(&la).then(a.steps[0].unitig.cmp(&b.steps[0].unitig))
    });
    Ok(layouts)
}

/// Extends the layout from `cur` through its free end `cur_end`.
///
/// `right` selects the chain side: `true` appends to the end, `false`
/// prepends. `cur_strand` is the orientation of `cur` in the contig.
fn extend(
    steps: &mut Vec<LayoutStep>,
    ends: &[Vec<Edge>],
    placed: &mut [bool],
    cur: usize,
    cur_end: usize,
    cur_strand: char,
    right: bool,
) {
    let mut cur = cur;
    let mut cur_end = cur_end;
    let mut cur_strand = cur_strand;
    loop {
        if is_repeat(ends, cur, cur_end) {
            break;
        }
        let Some(e) = best_edge(ends, cur, cur_end) else {
            break;
        };
        if placed[e.to] {
            break;
        }
        // Mutual-best junction: the target's best edge at the junction end
        // must point back to us (avoids joining through ambiguous branches).
        if best_edge(ends, e.to, e.to_end).map(|b| b.to) != Some(cur) {
            break;
        }
        placed[e.to] = true;
        let strand = if e.flip {
            flip_strand(cur_strand)
        } else {
            cur_strand
        };
        let mut step = LayoutStep {
            unitig: e.to,
            strand,
            q_start: 0,
            q_end: 0,
            overlap_len: 0,
        };
        if right {
            step.overlap_len = e.length;
            steps.push(step);
        } else {
            // The new leftmost step has no previous overlap; the former
            // first step now overlaps it by `e.length`.
            if let Some(first) = steps.first_mut() {
                first.overlap_len = e.length;
            }
            steps.insert(0, step);
        }
        cur = e.to;
        cur_end = 1 - e.to_end;
        cur_strand = strand;
    }
}

/// Pushes one directed edge (deduped per (to, to_end, flip) keeping the
/// longest length).
fn push_edge(
    ends: &mut [Vec<Edge>],
    from: usize,
    from_end: usize,
    to: usize,
    to_end: usize,
    flip: bool,
    length: usize,
) {
    ends[2 * from + from_end].push(Edge {
        from,
        from_end,
        to,
        to_end,
        flip,
        length,
    });
}

/// Deduplicates edges per target and sorts by (length desc, to, to_end,
/// flip) for a deterministic best/second-best.
fn dedup_and_sort(edges: &mut Vec<Edge>) {
    edges.sort_by(|a, b| {
        b.length
            .cmp(&a.length)
            .then(a.to.cmp(&b.to))
            .then(a.to_end.cmp(&b.to_end))
            .then(a.flip.cmp(&b.flip))
    });
    edges.dedup_by(|a, b| a.to == b.to && a.to_end == b.to_end && a.flip == b.flip);
}

/// Best extension edge at a unitig end.
fn best_edge(ends: &[Vec<Edge>], u: usize, end: usize) -> Option<Edge> {
    ends[2 * u + end].first().copied()
}

/// True when the unitig end has two near-equal best edges to different
/// targets (ambiguous junction; repeat evidence).
fn is_repeat(ends: &[Vec<Edge>], u: usize, end: usize) -> bool {
    let edges = &ends[2 * u + end];
    edges
        .first()
        .zip(edges.get(1))
        .is_some_and(|(best, second)| second.length * 10 >= best.length * 9 && second.to != best.to)
}

/// Opposite strand.
fn flip_strand(s: char) -> char {
    if s == '+' {
        '-'
    } else {
        '+'
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::nt::rev_comp;
    use crate::libs::olc::overlap::{find_overlaps, OverlapOptions};

    fn unitigs(names: &[&str], seqs: &[&str]) -> Vec<Unitig> {
        names
            .iter()
            .zip(seqs)
            .map(|(n, s)| Unitig {
                name: (*n).to_string(),
                seq: s.as_bytes().to_vec(),
            })
            .collect()
    }

    fn overlaps(us: &[Unitig], seed_k: usize, min_overlap: usize) -> Vec<Overlap> {
        find_overlaps(
            us,
            &OverlapOptions {
                seed_k,
                min_overlap,
            },
        )
        .unwrap()
    }

    /// Three unitigs chained by 8 bp suffix/prefix overlaps.
    #[test]
    fn linear_plus_chain() {
        let us = unitigs(
            &["u0", "u1", "u2"],
            &["AAAAAAAAACGTACGT", "ACGTACGTCCCCCCCC", "CCCCCCCCGGGGGGGG"],
        );
        let layouts = build_layouts(&us, &overlaps(&us, 5, 8)).unwrap();
        assert_eq!(layouts.len(), 1, "one linear chain");
        let l = &layouts[0];
        assert_eq!(l.steps.len(), 3);
        assert_eq!(l.steps[0].unitig, 0);
        assert_eq!(l.steps[0].strand, '+');
        assert_eq!((l.steps[0].q_start, l.steps[0].q_end), (0, 16));
        assert_eq!(l.steps[1].unitig, 1);
        assert_eq!(l.steps[1].strand, '+');
        assert_eq!((l.steps[1].q_start, l.steps[1].q_end), (8, 24));
        assert_eq!(l.steps[1].overlap_len, 8);
        assert_eq!(l.steps[2].unitig, 2);
        assert_eq!((l.steps[2].q_start, l.steps[2].q_end), (16, 32));
        assert_eq!(l.steps[2].overlap_len, 8);
    }

    /// A reverse-complement continuation: rc(u1) extends u0 to the right.
    #[test]
    fn linear_reverse_chain() {
        let u0 = "TTTTACGTAC";
        // rc(u1) = "ACGTAC" + "CCCC" so u1 = rc("ACGTACCCCC").
        let u1 = String::from_utf8(rev_comp(b"ACGTACCCCC").collect()).unwrap();
        let us = unitigs(&["u0", "u1"], &[u0, &u1]);
        let layouts = build_layouts(&us, &overlaps(&us, 5, 6)).unwrap();
        assert_eq!(layouts.len(), 1, "one reverse chain");
        let l = &layouts[0];
        assert_eq!(l.steps.len(), 2);
        assert_eq!(l.steps[0].unitig, 0);
        assert_eq!(l.steps[0].strand, '+');
        assert_eq!((l.steps[0].q_start, l.steps[0].q_end), (0, 10));
        assert_eq!(l.steps[1].unitig, 1);
        assert_eq!(l.steps[1].strand, '-');
        assert_eq!((l.steps[1].q_start, l.steps[1].q_end), (4, 14));
        assert_eq!(l.steps[1].overlap_len, 6);
    }

    /// Two equally good extensions stop at the branch (repeat junction).
    #[test]
    fn branch_stops_at_ambiguity() {
        let us = unitigs(
            &["u0", "u1", "u2"],
            &["AAAAACGTACGT", "ACGTACGTCCCC", "ACGTACGTGGGG"],
        );
        let layouts = build_layouts(&us, &overlaps(&us, 5, 8)).unwrap();
        assert_eq!(layouts.len(), 3, "no chain through the branch");
        assert!(layouts.iter().all(|l| l.steps.len() == 1));
    }

    /// Non-reciprocal junctions are not joined: u0's suffix matches u1's
    /// prefix, but u1's best edge at its prefix is a longer overlap to X.
    #[test]
    fn non_reciprocal_junction_not_joined() {
        let us = unitigs(
            &["u0", "u1", "u2"],
            &[
                "AAAAACGTACGT",     // suffix 8 = ACGTACGT
                "ACGTACGTACGTCCCC", // prefix 12 = ACGTACGTACGT
                "TTTTACGTACGTACGT", // suffix 12 = ACGTACGTACGT
            ],
        );
        let ovs = overlaps(&us, 5, 8);
        // u0->u1 (8 bp) and u2->u1 (12 bp) must both exist.
        assert!(ovs
            .iter()
            .any(|o| o.qid == 0 && o.tid == 1 && o.length == 8));
        assert!(ovs
            .iter()
            .any(|o| o.qid == 2 && o.tid == 1 && o.length == 12));
        let layouts = build_layouts(&us, &ovs).unwrap();
        // u2 (16 bp) joins u1; u0 (12 bp) cannot join u1 and stays alone.
        assert_eq!(layouts.len(), 2);
        let big = &layouts[0];
        assert_eq!(big.steps.len(), 2);
        assert_eq!(big.steps[0].unitig, 2);
        assert_eq!(big.steps[1].unitig, 1);
        assert_eq!(layouts[1].steps.len(), 1);
        assert_eq!(layouts[1].steps[0].unitig, 0);
    }

    /// Contain overlaps never create extension edges.
    #[test]
    fn contain_overlaps_do_not_chain() {
        let us = unitigs(
            &["long", "short"],
            &["AAAAGGTTAACCGGTTCCCC", "GGTTAACCGGTT"],
        );
        let ovs = overlaps(&us, 5, 10);
        assert!(ovs
            .iter()
            .any(|o| o.otype == super::super::overlap::OverlapType::Contain));
        let layouts = build_layouts(&us, &ovs).unwrap();
        assert!(layouts.iter().all(|l| l.steps.len() == 1));
    }

    /// The seed is the longest unitig and extends in both directions; the
    /// prepended (left) step's coordinates must be filled from its own
    /// length, not the seed's placeholder (regression: overflow panic).
    #[test]
    fn seed_extends_both_directions() {
        let us = unitigs(
            &["u0", "u1", "u2"],
            &[
                "AAAAACGTACGT",     // suffix 8 = ACGTACGT
                "ACGTACGTCCCCCCCC", // prefix 8 = ACGTACGT, suffix 8 = CCCCCCCC
                "CCCCCCCCGGGGGGGG", // prefix 8 = CCCCCCCC
            ],
        );
        let layouts = build_layouts(&us, &overlaps(&us, 5, 8)).unwrap();
        assert_eq!(layouts.len(), 1);
        let l = &layouts[0];
        assert_eq!(l.steps.len(), 3);
        assert_eq!(l.steps[0].unitig, 0);
        assert_eq!((l.steps[0].q_start, l.steps[0].q_end), (0, 12));
        assert_eq!(l.steps[0].overlap_len, 0);
        assert_eq!(l.steps[1].unitig, 1);
        assert_eq!((l.steps[1].q_start, l.steps[1].q_end), (4, 20));
        assert_eq!(l.steps[1].overlap_len, 8);
        assert_eq!(l.steps[2].unitig, 2);
        assert_eq!((l.steps[2].q_start, l.steps[2].q_end), (12, 28));
        assert_eq!(l.steps[2].overlap_len, 8);
    }

    /// An overlap longer than the previous step (malformed user PAF) is a
    /// friendly error, not a panic (zero-panic policy).
    #[test]
    fn inconsistent_overlap_is_error() {
        let us = unitigs(&["u0", "u1"], &["AAAACCCC", "CCCCGGGG"]);
        let ovs = vec![Overlap {
            qid: 0,
            tid: 1,
            strand: '+',
            q_start: 0,
            q_end: 8,
            t_start: 0,
            t_end: 8,
            length: 20,
            otype: crate::libs::olc::overlap::OverlapType::Dovetail,
        }];
        let err = build_layouts(&us, &ovs).unwrap_err();
        assert!(err.to_string().contains("exceeds"), "{err}");
    }
}
