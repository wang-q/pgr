//! Window derivation from FasBlock reference entries.
//!
//! [`derive_windows_from_blocks`] scans reference ranges across all input
//! block sets, unions overlapping intervals (expanded by `radius`), and keeps
//! windows satisfying the `min_width` and coverage requirements.

use super::{find_ref_entry, FasMultizConfig, Window};
use crate::libs::ds::merge_intervals;
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
                if range.start() > range.end() {
                    // Inverted (malformed) reference range, e.g. `>ref.chr(+):100-1`.
                    // Deriving an interval from it would make `width = e - s`
                    // below underflow (debug panic, release wrap). Skip it.
                    continue;
                }
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

    // No coverage filter is needed: every window is derived from (at least
    // one) input's reference interval expanded by `radius`, so it is always
    // covered by that input. A per-input DupeTree filter here would be
    // redundant for normal windows and would silently drop single-base
    // reference blocks (whose zero-width interval `DupeTree::add` ignores).
    windows
}
