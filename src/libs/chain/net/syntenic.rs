//! Synteny classification: depth-tracking interval tree + net walker.
//!
//! `DupeTree` accumulates signed intervals (added by fills, subtracted by
//! nested gaps) and flattens them into non-overlapping `Segment`s of constant
//! depth, so a fill's query-overlap with duplications can be queried.
//!
//! `classify_syntenic` walks a Net's gap-fill tree and assigns each fill a
//! synteny class (`top`/`syn`/`inv`/`nonSyn`) plus qOver/qFar/qDup statistics.

use crate::libs::chain::net::writer::range_intersection;
use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

use super::types::{Chrom, Fill, Gap};

/// A non-overlapping run of constant duplication depth.
pub struct Segment {
    pub start: u64,
    pub end: u64,
    pub depth: i32,
}

/// Interval tree tracking query-side duplication depth for synteny classification.
pub struct DupeTree {
    intervals: Vec<(u64, u64, i32)>,
    segments: Vec<Segment>,
}

impl Default for DupeTree {
    fn default() -> Self {
        Self::new()
    }
}

impl DupeTree {
    /// Creates an empty DupeTree.
    pub fn new() -> Self {
        Self {
            intervals: Vec::new(),
            segments: Vec::new(),
        }
    }

    /// Records a +1 depth contribution over `[start, end)`.
    pub fn add(&mut self, start: u64, end: u64) {
        if start < end {
            self.intervals.push((start, end, 1));
        }
    }

    /// Records a -1 depth contribution over `[start, end)`.
    pub fn subtract(&mut self, start: u64, end: u64) {
        if start < end {
            self.intervals.push((start, end, -1));
        }
    }

    /// Flattens recorded intervals into non-overlapping constant-depth segments.
    pub fn build(&mut self) {
        if self.intervals.is_empty() {
            return;
        }

        let mut events = Vec::new();
        for (s, e, d) in &self.intervals {
            events.push((*s, *d));
            events.push((*e, -*d));
        }
        // Sort by position, then by delta so that end events (negative delta)
        // are processed before start events (positive delta) at the same coordinate.
        events.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));

        let mut current_depth = 0;
        let mut segments = Vec::new();

        for i in 0..events.len() - 1 {
            let (pos, delta) = events[i];
            current_depth += delta;

            let next_pos = events[i + 1].0;
            if next_pos > pos {
                segments.push(Segment {
                    start: pos,
                    end: next_pos,
                    depth: current_depth,
                });
            }
        }

        self.segments = segments;
    }

    /// Returns total bases covered by segments with `depth >= threshold`.
    pub fn count_over(&self, start: u64, end: u64, threshold: i32) -> u64 {
        if self.segments.is_empty() {
            return 0;
        }

        // Binary search for first segment ending after start
        let idx = self.segments.partition_point(|seg| seg.end <= start);

        let mut count = 0;
        for seg in &self.segments[idx..] {
            if seg.start >= end {
                break;
            }

            if seg.depth >= threshold {
                let s = seg.start.max(start);
                let e = seg.end.min(end);
                if s < e {
                    count += e - s;
                }
            }
        }
        count
    }
}

/// Classify syntenic relationship for all fills in `chroms`.
///
/// For each query chromosome referenced by fills, builds a `DupeTree`
/// tracking signed depth contributions (fills add, nested gaps subtract),
/// then walks each fill assigning:
/// * `class` — `top` (root fill), `syn`/`inv` (same query chrom, same/opposite strand),
///   `nonSyn` (different query chrom from parent fill)
/// * `q_dup` — bases of the fill's query range with duplication depth >= 2
/// * `q_over` — bases of overlap with the parent gap's query range
/// * `q_far` — distance from the fill's query range to the parent gap's query range
///   (0 when they overlap)
pub fn classify_syntenic(chroms: &[Chrom]) {
    // Build DupeTrees for all query chromosomes
    let mut q_chrom_map: HashMap<String, DupeTree> = HashMap::new();
    for chrom in chroms {
        r_calc_dupes(chrom, &mut q_chrom_map);
    }

    // Flatten intervals into constant-depth segments
    for dt in q_chrom_map.values_mut() {
        dt.build();
    }

    // Classify every fill in the net tree
    for chrom in chroms {
        r_net_syn(chrom, &q_chrom_map);
    }
}

fn r_calc_dupes(chrom: &Chrom, map: &mut HashMap<String, DupeTree>) {
    r_calc_dupes_gap(&chrom.root, map);
}

fn r_calc_dupes_gap(gap: &Rc<RefCell<Gap>>, map: &mut HashMap<String, DupeTree>) {
    let g = gap.borrow();
    for fill in &g.fills {
        r_calc_dupes_fill(fill, map);
    }
}

fn r_calc_dupes_fill(fill: &Rc<RefCell<Fill>>, map: &mut HashMap<String, DupeTree>) {
    let f = fill.borrow();
    let q_name = &f.o_chrom;
    let start = f.o_start;
    let end = f.o_end;

    if !q_name.is_empty() {
        let dt = map.entry(q_name.clone()).or_default();
        dt.add(start, end);
    }

    // Recursively process gaps inside fill
    for gap in &f.gaps {
        let g = gap.borrow();
        // Gap inside Fill shares query chrom with Fill
        // But Gap subtracts coverage
        let q_name = &f.o_chrom;
        let start = g.o_start;
        let end = g.o_end;

        if !q_name.is_empty() {
            let dt = map.entry(q_name.clone()).or_default();
            dt.subtract(start, end);
        }

        // Recurse into fills inside gap
        r_calc_dupes_gap(gap, map);
    }
}

fn r_net_syn(chrom: &Chrom, map: &HashMap<String, DupeTree>) {
    r_net_syn_gap(&chrom.root, map, None);
}

fn r_net_syn_gap(
    gap: &Rc<RefCell<Gap>>,
    map: &HashMap<String, DupeTree>,
    parent_fill: Option<&Rc<RefCell<Fill>>>,
) {
    let g = gap.borrow();
    for fill in &g.fills {
        r_net_syn_fill(fill, map, parent_fill, Some(gap));
    }
}

fn r_net_syn_fill(
    fill: &Rc<RefCell<Fill>>,
    map: &HashMap<String, DupeTree>,
    parent: Option<&Rc<RefCell<Fill>>>,
    parent_gap: Option<&Rc<RefCell<Gap>>>,
) {
    // Need to borrow_mut to update fields
    // But we also need to pass `fill` (Rc) to children.
    // So we borrow mut, update, drop borrow, then recurse.

    let (q_name, start, end, strand) = {
        let f = fill.borrow();
        (f.o_chrom.clone(), f.o_start, f.o_end, f.o_strand)
    };

    let q_dup = if let Some(dt) = map.get(&q_name) {
        Some(dt.count_over(start, end, 2))
    } else {
        Some(0)
    };

    let mut q_over = None;
    let mut q_far = None;

    let type_str = match parent {
        None => "top".to_string(),
        Some(p_rc) => {
            let p = p_rc.borrow();
            if q_name != p.o_chrom {
                "nonSyn".to_string()
            } else {
                // Calculate qOver/qFar relative to parent GAP
                if let Some(pg_rc) = parent_gap {
                    let pg = pg_rc.borrow();
                    // Check overlap with GAP query range
                    let g_start = pg.o_start;
                    let g_end = pg.o_end;

                    let intersection = range_intersection(start, end, g_start, g_end);
                    q_over = Some(intersection);

                    if intersection == 0 {
                        // Calculate distance
                        let d1 = start.saturating_sub(g_end);
                        let d2 = g_start.saturating_sub(end);
                        q_far = Some((d1 + d2) as i64);
                    } else {
                        q_far = Some(0);
                    }
                } else {
                    // Should not happen for non-top fills
                    q_over = Some(0);
                    q_far = Some(0);
                }

                if p.o_strand == strand {
                    "syn".to_string()
                } else {
                    "inv".to_string()
                }
            }
        }
    };

    {
        let mut f = fill.borrow_mut();
        f.class = type_str;
        f.q_dup = q_dup;
        f.q_over = q_over;
        f.q_far = q_far;
    }

    // Recurse
    // Children of fill are in `f.gaps`
    // We need to access `f.gaps` without holding mutable borrow on `f`
    let gaps = fill.borrow().gaps.clone();
    for gap in gaps {
        r_net_syn_gap(&gap, map, Some(fill));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dupe_tree_adjacent_intervals() {
        let mut dt = DupeTree::new();
        dt.add(0, 10);
        dt.add(10, 20);
        dt.build();

        // Each adjacent interval contributes depth 1 on its own range.
        assert_eq!(dt.count_over(0, 20, 1), 20);
        assert_eq!(dt.count_over(0, 20, 2), 0);
    }

    #[test]
    fn test_dupe_tree_overlapping_intervals() {
        let mut dt = DupeTree::new();
        dt.add(0, 15);
        dt.add(10, 25);
        dt.build();

        assert_eq!(dt.count_over(0, 25, 1), 25);
        assert_eq!(dt.count_over(0, 25, 2), 5);
        assert_eq!(dt.count_over(0, 10, 2), 0);
        assert_eq!(dt.count_over(20, 25, 2), 0);
    }
}
