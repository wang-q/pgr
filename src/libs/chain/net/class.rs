//! Net class statistics: aggregate fill/gap bases by class label.
//!
//! Walks the net tree recursively, tallying per-class count and bases for
//! fills and the implicit gap background (gap size minus enclosed fills).

use super::types::{Fill, Gap};
use std::cell::{Ref, RefCell};
use std::collections::HashMap;
use std::rc::Rc;

/// Per-class aggregate counters (count and bases covered).
#[derive(Default)]
pub struct Stats {
    pub count: u64,
    pub bases: u64,
}

/// Recursively collect stats for a gap and its nested fills.
pub fn collect_stats_gap(gap: &Rc<RefCell<Gap>>, stats: &mut HashMap<String, Stats>) {
    let gap_ref: Ref<Gap> = gap.borrow();
    let size = gap_ref.end - gap_ref.start;

    // Gap itself is a "gap" class if we want to count it?
    // Or do we only count fills?
    // UCSC netClass counts "gap" as well.
    // But Gaps contain Fills.
    // The "gap" bases are (size - sum(fills.size)).

    // Actually, usually we count the explicit objects.
    // A Gap object represents a gap in the alignment.
    // But in the net structure, Gap is a container.
    // The "unfilled" part of the Gap is the actual gap.

    let mut fill_bases = 0;
    for fill in &gap_ref.fills {
        let fill_ref: Ref<Fill> = fill.borrow();
        fill_bases += fill_ref.end - fill_ref.start;

        // Count the fill
        let class = if fill_ref.class.is_empty() {
            "unknown".to_string()
        } else {
            fill_ref.class.clone()
        };

        let entry = stats.entry(class).or_default();
        entry.count += 1;
        entry.bases += fill_ref.end - fill_ref.start;

        // Recurse
        collect_stats_fill(fill, stats);
    }

    // The remaining part is gap
    let gap_bases = size.saturating_sub(fill_bases);
    if gap_bases > 0 {
        let entry = stats.entry("gap".to_string()).or_default();
        entry.count += 1; // This is tricky. Is it 1 gap? Or multiple fragments?
                          // In this recursive structure, the "gap" is implicitly the background.
                          // We can just add the bases.
                          // Count is hard to define for the implicit background gap.
                          // Let's just track bases for gap.
        entry.bases += gap_bases;
    }
}

/// Recursively collect stats for a fill's nested gaps.
pub fn collect_stats_fill(fill: &Rc<RefCell<Fill>>, stats: &mut HashMap<String, Stats>) {
    let fill_ref: Ref<Fill> = fill.borrow();

    // Fill contains Gaps.
    // The Fill itself covers fill_ref.end - fill_ref.start.
    // This was already added to the stats in the parent.
    // But we need to recurse into its gaps.

    for gap in &fill_ref.gaps {
        collect_stats_gap(gap, stats);
    }
}
