//! Window derivation from FasBlock reference entries.
//!
//! [`derive_windows_from_blocks`] scans reference ranges across all input
//! block sets, unions overlapping intervals (expanded by `radius`), and keeps
//! windows satisfying the `min_width` and per-mode coverage requirements.

use super::{find_ref_entry, ref_overlaps_window, FasMultizConfig, FasMultizMode, Window};
use crate::libs::fmt::fas::FasBlock;
use std::collections::BTreeMap;

pub(super) fn derive_windows_from_blocks(
    ref_name: &str,
    blocks_per_input: &[Vec<FasBlock>],
    cfg: &FasMultizConfig,
) -> Vec<Window> {
    use std::cmp::max;

    let mut per_chr: BTreeMap<String, Vec<(u64, u64)>> = BTreeMap::new();

    for group in blocks_per_input {
        for block in group {
            if let Some(entry) = find_ref_entry(block, ref_name) {
                let range = entry.range();
                let chr = range.chr().to_string();
                let start = *range.start() as u64;
                let end = *range.end() as u64;
                let s = start.saturating_sub(cfg.radius as u64);
                let e = end + cfg.radius as u64;
                per_chr.entry(chr).or_default().push((s, e));
            }
        }
    }

    let mut windows = Vec::new();

    for (chr, mut intervals) in per_chr {
        if intervals.is_empty() {
            continue;
        }
        intervals.sort_by_key(|(s, _)| *s);

        let mut current = intervals[0];
        for &(s, e) in &intervals[1..] {
            if s <= current.1 {
                current.1 = max(current.1, e);
            } else {
                if current.1 > current.0 {
                    let width = current.1 - current.0;
                    if width >= cfg.min_width as u64 {
                        windows.push(Window {
                            chr: chr.clone(),
                            start: current.0,
                            end: current.1,
                        });
                    }
                }
                current = (s, e);
            }
        }

        if current.1 > current.0 {
            let width = current.1 - current.0;
            if width >= cfg.min_width as u64 {
                windows.push(Window {
                    chr,
                    start: current.0,
                    end: current.1,
                });
            }
        }
    }

    if windows.is_empty() {
        return windows;
    }

    let total_inputs = blocks_per_input.len();
    let required_inputs = match cfg.mode {
        FasMultizMode::Core => total_inputs,
        FasMultizMode::Union => 1,
    };

    let mut filtered = Vec::new();

    for window in windows {
        let mut covered = 0usize;
        for group in blocks_per_input {
            let has_overlap = group.iter().any(|block| {
                find_ref_entry(block, ref_name)
                    .map(|entry| ref_overlaps_window(entry, &window))
                    .unwrap_or(false)
            });
            if has_overlap {
                covered += 1;
            }
            if covered >= required_inputs {
                break;
            }
        }

        if covered >= required_inputs {
            filtered.push(window);
        }
    }

    filtered
}
