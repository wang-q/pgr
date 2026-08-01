//! Window derivation from FasBlock reference entries.
//!
//! [`derive_windows_from_blocks`] scans reference ranges across all input
//! block sets, unions overlapping intervals (expanded by `radius`), and keeps
//! windows satisfying the `min_width` and per-mode coverage requirements.

use super::{find_ref_entry, FasMultizConfig, FasMultizMode, Window};
use crate::libs::ds::{merge_intervals, DupeTree};
use crate::libs::fmt::fas::FasBlock;
use std::collections::BTreeMap;

pub(super) fn derive_windows_from_blocks(
    ref_name: &str,
    blocks_per_input: &[Vec<FasBlock>],
    cfg: &FasMultizConfig,
) -> Vec<Window> {
    let mut per_chr: BTreeMap<String, Vec<(u64, u64)>> = BTreeMap::new();

    for group in blocks_per_input {
        for block in group {
            if let Some(entry) = find_ref_entry(block, ref_name) {
                let range = entry.range();
                let chr = range.chr().to_string();
                let start = *range.start() as u64;
                let end = *range.end() as u64;
                let s = start.saturating_sub(cfg.radius as u64);
                let e = end.saturating_add(cfg.radius as u64);
                per_chr.entry(chr).or_default().push((s, e));
            }
        }
    }

    let mut windows = Vec::new();

    for (chr, mut intervals) in per_chr {
        if intervals.is_empty() {
            continue;
        }
        for (s, e) in merge_intervals(&mut intervals) {
            let width = e - s;
            if width >= cfg.min_width as u64 {
                windows.push(Window {
                    chr: chr.clone(),
                    start: s,
                    end: e,
                });
            }
        }
    }

    if windows.is_empty() {
        return windows;
    }

    let required_inputs = match cfg.mode {
        FasMultizMode::Core => blocks_per_input.len() as i32,
        FasMultizMode::Union => 1,
    };

    // Per-chromosome DupeTree: each input contributes at most 1 depth over its
    // (merged) reference intervals, so `count_over(window, required) > 0`
    // means the window overlaps at least `required` distinct inputs.
    let mut cov_trees: BTreeMap<String, DupeTree> = BTreeMap::new();
    for group in blocks_per_input {
        let mut by_chr: BTreeMap<String, Vec<(u64, u64)>> = BTreeMap::new();
        for block in group {
            if let Some(entry) = find_ref_entry(block, ref_name) {
                let range = entry.range();
                by_chr
                    .entry(range.chr().to_string())
                    .or_default()
                    .push((*range.start() as u64, *range.end() as u64));
            }
        }
        for (chr, mut intervals) in by_chr {
            let tree = cov_trees.entry(chr).or_default();
            for (s, e) in merge_intervals(&mut intervals) {
                tree.add(s, e);
            }
        }
    }
    for tree in cov_trees.values_mut() {
        tree.build();
    }

    let mut filtered = Vec::new();
    for window in windows {
        let covered = cov_trees.get(&window.chr).map_or(0, |tree| {
            tree.count_over(window.start, window.end, required_inputs)
        });
        if covered > 0 {
            filtered.push(window);
        }
    }

    filtered
}
