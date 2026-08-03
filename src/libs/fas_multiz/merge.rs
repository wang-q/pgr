//! FasBlock merging: pairwise DP merge and per-window block assembly.

use super::banded_align::banded_align_refs;
use super::{find_ref_entry, ref_overlaps_window, FasMultizConfig, Window};
use crate::libs::chain::sub_matrix::SubMatrix;
use crate::libs::ds::best_crossover;
use crate::libs::fmt::fas::{FasBlock, FasEntry};
use std::collections::{BTreeMap, BTreeSet};

/// Sorted species name set of a block.
fn block_species(block: &FasBlock) -> BTreeSet<String> {
    block.names.iter().cloned().collect()
}

/// Deterministic content key for ordering blocks: reference range plus the
/// sorted species list. Depends only on block contents, never on input order.
fn content_key(block: &FasBlock, ref_name: &str) -> (String, u64, u64, String) {
    match find_ref_entry(block, ref_name) {
        Some(entry) => {
            let range = entry.range();
            (
                range.chr().to_string(),
                *range.start() as u64,
                *range.end() as u64,
                block.names.join(","),
            )
        }
        None => (String::new(), 0, 0, block.names.join(",")),
    }
}

/// Pick the candidate with the largest species overlap with `union`, breaking
/// ties by species count and then by content key.
fn best_block(
    blocks: &[&FasBlock],
    species: &[BTreeSet<String>],
    union: &BTreeSet<String>,
    candidates: &[usize],
    ref_name: &str,
) -> usize {
    let priority = |i: usize| {
        (
            species[i].intersection(union).count(),
            species[i].len(),
            content_key(blocks[i], ref_name),
        )
    };
    let mut best = candidates[0];
    for &idx in &candidates[1..] {
        if priority(idx) > priority(best) {
            best = idx;
        }
    }
    best
}

/// Deterministic content-based merge order for progressive DP merging.
///
/// Greedy agglomeration: start with the block carrying the most species, then
/// repeatedly attach the remaining block with the largest species overlap to
/// the accumulated union. All tie-breaks use block contents, so the order is
/// independent of the input file order.
fn merge_order(blocks: &[&FasBlock], ref_name: &str) -> Vec<usize> {
    let species: Vec<BTreeSet<String>> = blocks.iter().map(|b| block_species(b)).collect();
    let mut remaining: Vec<usize> = (0..blocks.len()).collect();
    let mut union: BTreeSet<String> = BTreeSet::new();
    let mut order = Vec::with_capacity(blocks.len());
    while !remaining.is_empty() {
        let idx = best_block(blocks, &species, &union, &remaining, ref_name);
        union.extend(species[idx].iter().cloned());
        remaining.retain(|&i| i != idx);
        order.push(idx);
    }
    order
}

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

fn find_species_entry<'a>(block: &'a FasBlock, name: &str) -> Option<&'a FasEntry> {
    block
        .entries
        .iter()
        .zip(block.names.iter())
        .find_map(|(entry, n)| if n == name { Some(entry) } else { None })
}

/// Try merging two blocks whose reference sequences differ beyond gap
/// placement: align the two column profiles, cut the overlap at the best
/// crossover point, and splice the left part from `blocks[0]` with the right
/// part from `blocks[1]` (mirrors the `best_crossover` idea from UCSC).
fn merge_conflicting_refs(
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

    // A shared non-reference species is required to score the two profiles.
    let shared = blocks[0]
        .names
        .iter()
        .find(|n| *n != ref_name && blocks[1].names.iter().any(|m| m == *n));
    let Some(shared_name) = shared else {
        return Ok(None);
    };
    let x_a = find_species_entry(blocks[0], shared_name);
    let x_b = find_species_entry(blocks[1], shared_name);
    let (Some(x_a), Some(x_b)) = (x_a, x_b) else {
        return Ok(None);
    };

    let (map_a, map_b) = match banded_align_refs(blocks, ref_name, cfg) {
        Some(v) => v,
        None => return Ok(None),
    };
    let out_len = map_a.len();

    // Hardcoded multiz scoring: HOX70 (= hoxd55) substitution matrix.
    let submat = SubMatrix::hoxd55();

    let col = |map: &[Option<usize>], seq: &[u8], i: usize| -> u8 {
        map[i].and_then(|idx| seq.get(idx).copied()).unwrap_or(b'-')
    };

    let mut l_t = Vec::with_capacity(out_len);
    let mut l_q = Vec::with_capacity(out_len);
    let mut r_t = Vec::with_capacity(out_len);
    let mut r_q = Vec::with_capacity(out_len);
    for i in 0..out_len {
        l_t.push(col(&map_a, ref_a.seq(), i));
        l_q.push(col(&map_a, x_a.seq(), i));
        r_t.push(col(&map_b, ref_b.seq(), i));
        r_q.push(col(&map_b, x_b.seq(), i));
    }

    let (cut, _) = best_crossover(&l_t, &l_q, &r_t, &r_q, |a, b| {
        submat.get_score(a as char, b as char) as f64
    });

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

    let mut entries = Vec::new();
    let mut names = Vec::new();
    let mut headers = Vec::new();

    for name in species {
        let group = species_map.get(&name).unwrap();

        let mut seq = Vec::with_capacity(out_len);
        for pos in 0..out_len {
            let (map, entry) = if pos < cut {
                (&map_a, group[0])
            } else {
                (&map_b, group[1])
            };
            let base = match entry {
                Some(e) => col(map, e.seq(), pos),
                None => b'-',
            };
            seq.push(base);
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
        return merge_conflicting_refs(ref_name, blocks, cfg);
    }

    let (map_a, map_b) = match banded_align_refs(blocks, ref_name, cfg) {
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

            // The merged reference preserves the first block's reference
            // sequence: filling ref_a gaps with ref_b bases would inflate the
            // ungapped length and break the `ungapped_equal` invariant that
            // downstream merges rely on. Non-reference species keep the
            // two-block fallback.
            if name != ref_name && chosen.is_none() {
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

    let order = merge_order(blocks, ref_name);
    let ordered: Vec<&FasBlock> = order.iter().map(|&i| blocks[i]).collect();

    let mut acc = match merge_two_blocks_with_dp(ref_name, [ordered[0], ordered[1]], cfg)? {
        Some(v) => v,
        None => return Ok(None),
    };

    if ordered.len() == 2 {
        return Ok(Some(acc));
    }

    for &block in &ordered[2..] {
        if let Some(next) = merge_two_blocks_with_dp(ref_name, [&acc, block], cfg)? {
            acc = next;
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
        if let Some(block) = candidate {
            blocks.push(block);
        } else {
            // Inputs without a block in this window are simply skipped.
        }
    }

    if blocks.is_empty() {
        return Ok(None);
    }

    // Deterministic order for the non-DP fallback below (first entry per
    // species must not depend on input file order).
    blocks.sort_by_key(|b| content_key(b, ref_name));

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
