//! FasBlock merging: pairwise DP merge and per-window block assembly.

use super::banded_align::banded_align_refs;
use super::{find_ref_entry, ref_overlaps_window, FasMultizConfig, FasMultizMode, Window};
use crate::libs::fmt::fas::{FasBlock, FasEntry};
use std::collections::BTreeMap;

fn entry_seq_equal(a: &FasEntry, b: &FasEntry) -> bool {
    a.seq() == b.seq()
}

fn ungapped_equal(a: &FasEntry, b: &FasEntry) -> bool {
    let sa = a.seq();
    let sb = b.seq();
    let ua: Vec<u8> = sa.iter().copied().filter(|c| *c != b'-').collect();
    let ub: Vec<u8> = sb.iter().copied().filter(|c| *c != b'-').collect();
    ua == ub
}

fn merge_two_blocks_with_dp(
    ref_name: &str,
    blocks: [&FasBlock; 2],
    cfg: &FasMultizConfig,
) -> anyhow::Result<Option<FasBlock>> {
    let ref_a = match find_ref_entry(blocks[0], ref_name) {
        Some(v) => v,
        None => return Ok(None),
    };
    let ref_b = match find_ref_entry(blocks[1], ref_name) {
        Some(v) => v,
        None => return Ok(None),
    };

    if !ungapped_equal(ref_a, ref_b) {
        return Ok(None);
    }

    let (map_a, map_b) = match banded_align_refs(blocks, ref_name, cfg)? {
        Some(v) => v,
        None => return Ok(None),
    };

    let ref_range = ref_a.range().clone();

    let mut species_map: BTreeMap<String, [Option<&FasEntry>; 2]> = BTreeMap::new();

    for (idx, block) in blocks.iter().enumerate() {
        for (entry, name) in block.entries.iter().zip(block.names.iter()) {
            let v = species_map.entry(name.clone()).or_insert([None, None]);
            v[idx] = Some(entry);
        }
    }

    let mut species: Vec<String> = species_map.keys().cloned().collect();
    species.sort();
    species.sort_by_key(|n| if n == ref_name { 0 } else { 1 });

    let out_len = map_a.len();

    let mut entries = Vec::new();
    let mut names = Vec::new();
    let mut headers = Vec::new();

    for name in species {
        let group = species_map.get(&name).unwrap();

        if matches!(cfg.mode, FasMultizMode::Core) && (group[0].is_none() || group[1].is_none()) {
            continue;
        }

        let mut seq = Vec::with_capacity(out_len);

        for pos in 0..out_len {
            let mut chosen: Option<u8> = None;

            if let Some(entry) = group[0] {
                if let Some(idx) = map_a[pos] {
                    if idx < entry.seq().len() {
                        chosen = Some(entry.seq()[idx]);
                    }
                }
            }

            if chosen.is_none() {
                if let Some(entry) = group[1] {
                    if let Some(idx) = map_b[pos] {
                        if idx < entry.seq().len() {
                            chosen = Some(entry.seq()[idx]);
                        }
                    }
                }
            }

            seq.push(chosen.unwrap_or(b'-'));
        }

        let range = if name == ref_name {
            ref_range.clone()
        } else {
            let chosen = if group[0].is_some() {
                group[0]
            } else {
                group[1]
            }
            .unwrap();
            chosen.range().clone()
        };

        let entry = FasEntry::from(&range, &seq);
        let header = format!("{}", range);

        entries.push(entry);
        names.push(name.clone());
        headers.push(header);
    }

    if entries.is_empty() {
        Ok(None)
    } else {
        Ok(Some(FasBlock {
            entries,
            names,
            headers,
        }))
    }
}

fn merge_blocks_with_dp(
    ref_name: &str,
    blocks: &[&FasBlock],
    cfg: &FasMultizConfig,
) -> anyhow::Result<Option<FasBlock>> {
    if blocks.len() < 2 {
        return Ok(None);
    }

    let mut acc = match merge_two_blocks_with_dp(ref_name, [blocks[0], blocks[1]], cfg)? {
        Some(v) => v,
        None => return Ok(None),
    };

    if blocks.len() == 2 {
        return Ok(Some(acc));
    }

    match cfg.mode {
        FasMultizMode::Core => {
            for &block in &blocks[2..] {
                acc = match merge_two_blocks_with_dp(ref_name, [&acc, block], cfg)? {
                    Some(v) => v,
                    None => return Ok(None),
                };
            }
        }
        FasMultizMode::Union => {
            for &block in &blocks[2..] {
                if let Some(next) = merge_two_blocks_with_dp(ref_name, [&acc, block], cfg)? {
                    acc = next;
                }
            }
        }
    }

    Ok(Some(acc))
}

pub fn merge_window(
    ref_name: &str,
    window: &Window,
    blocks_per_input: &[Vec<FasBlock>],
    cfg: &FasMultizConfig,
) -> anyhow::Result<Option<FasBlock>> {
    if blocks_per_input.is_empty() {
        return Ok(None);
    }

    let mut blocks = Vec::new();
    for group in blocks_per_input {
        let candidate = group
            .iter()
            .find(|block| match find_ref_entry(block, ref_name) {
                Some(entry) => ref_overlaps_window(entry, window),
                None => false,
            });
        match candidate {
            Some(block) => blocks.push(block),
            None => {
                if matches!(cfg.mode, FasMultizMode::Core) {
                    return Ok(None);
                }
            }
        }
    }

    if blocks.is_empty() {
        return Ok(None);
    }

    if blocks.len() >= 2 {
        if let Some(block) = merge_blocks_with_dp(ref_name, &blocks, cfg)? {
            return Ok(Some(block));
        }
    }

    let template = blocks[0];
    let ref_entry = match find_ref_entry(template, ref_name) {
        Some(v) => v,
        None => return Ok(None),
    };

    for block in &blocks[1..] {
        let other_ref = match find_ref_entry(block, ref_name) {
            Some(v) => v,
            None => return Ok(None),
        };
        if !entry_seq_equal(ref_entry, other_ref) {
            return Ok(None);
        }
    }

    let ref_range = ref_entry.range().clone();

    let n = blocks.len();
    let mut species_map: BTreeMap<String, Vec<Option<&FasEntry>>> = BTreeMap::new();

    for (i, block) in blocks.iter().enumerate() {
        for (entry, name) in block.entries.iter().zip(block.names.iter()) {
            let v = species_map
                .entry(name.clone())
                .or_insert_with(|| vec![None; n]);
            v[i] = Some(entry);
        }
    }

    let mut species: Vec<String> = species_map.keys().cloned().collect();
    species.sort();
    species.sort_by_key(|n| if n == ref_name { 0 } else { 1 });

    let mut entries = Vec::new();
    let mut names = Vec::new();
    let mut headers = Vec::new();

    for name in species {
        let group = species_map.get(&name).unwrap();

        if matches!(cfg.mode, FasMultizMode::Core) && group.iter().any(|e| e.is_none()) {
            continue;
        }

        let chosen = match group.iter().flatten().next() {
            Some(e) => e,
            None => continue,
        };

        let range = if name == ref_name {
            ref_range.clone()
        } else {
            chosen.range().clone()
        };

        let seq = chosen.seq();
        let entry = FasEntry::from(&range, seq);
        let header = format!("{}", range);

        entries.push(entry);
        names.push(name.clone());
        headers.push(header);
    }

    if entries.is_empty() {
        Ok(None)
    } else {
        Ok(Some(FasBlock {
            entries,
            names,
            headers,
        }))
    }
}
