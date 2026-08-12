//! FasBlock merging: pairwise DP merge and per-window block assembly.

use super::banded_align::banded_align_refs;
use super::{find_ref_entry, ref_overlaps_window, FasMultizConfig, Window};
use crate::libs::chain::sub_matrix::SubMatrix;
use crate::libs::ds::{best_crossover, Range};
use crate::libs::fmt::fas::{FasBlock, FasEntry};
use std::collections::BTreeMap;

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
        // A species present in only one block must be carried across the whole
        // merged output through that block's map, otherwise the half of the
        // sequence on the other side of `cut` collapses to gaps and the species
        // silently loses data. Shared species (and the reference) are spliced at
        // the crossover so the left part tracks block A and the right part
        // block B, matching the spliced reference.
        let splice = group[0].is_some() && group[1].is_some();
        for pos in 0..out_len {
            let base = if splice {
                let (map, entry) = if pos < cut {
                    (&map_a, group[0])
                } else {
                    (&map_b, group[1])
                };
                match entry {
                    Some(e) => col(map, e.seq(), pos),
                    None => b'-',
                }
            } else if group[0].is_some() {
                col(&map_a, group[0].unwrap().seq(), pos)
            } else {
                col(&map_b, group[1].unwrap().seq(), pos)
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
            // sequence; non-reference species keep the two-block fallback.
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

/// Reference range `(start, end)` of a block's reference entry.
fn ref_range(block: &FasBlock, ref_name: &str) -> Option<(i32, i32)> {
    let entry = find_ref_entry(block, ref_name)?;
    Some((*entry.range().start(), *entry.range().end()))
}

/// Reference size of a block (reference positions covered), 0 when the block
/// has no reference entry (e.g. an insertion-only part).
fn ref_size(block: &FasBlock, ref_name: &str) -> i32 {
    ref_range(block, ref_name).map_or(0, |(s, e)| e - s + 1)
}

/// Column index of the reference base at 1-based position `pos` (canonical
/// multiz `mafPos2Col`).
fn ref_pos_to_col(entry: &FasEntry, pos: i32) -> Option<usize> {
    let mut p = *entry.range().start() - 1;
    for (col, &base) in entry.seq().iter().enumerate() {
        if base != b'-' {
            p += 1;
            if p == pos {
                return Some(col);
            }
        }
    }
    None
}

/// Slice a block to columns `[cbeg, cend]`: recompute every entry's range
/// from the bases kept, drop species left with no bases, and remove all-gap
/// columns (canonical `make_part_ali_col` + `mafColDashRm`).
fn slice_block_cols(
    block: &FasBlock,
    ref_name: &str,
    cbeg: usize,
    cend: usize,
    ref_beg: i32,
    ref_end: i32,
) -> Option<FasBlock> {
    if cend < cbeg {
        return None;
    }
    let mut kept: Vec<(String, String, Range, Vec<u8>)> = Vec::new();
    for (i, (name, entry)) in block.names.iter().zip(block.entries.iter()).enumerate() {
        let slice = &entry.seq()[cbeg..=cend];
        let bases_before = entry.seq()[..cbeg].iter().filter(|&&b| b != b'-').count() as i32;
        let size = slice.iter().filter(|&&b| b != b'-').count();
        if size == 0 {
            continue;
        }
        let range = if name == ref_name {
            Range::from_full(
                entry.range().name().as_str(),
                entry.range().chr().as_str(),
                entry.range().strand().as_str(),
                ref_beg,
                ref_end,
            )
        } else {
            let start = *entry.range().start() + bases_before;
            Range::from_full(
                entry.range().name().as_str(),
                entry.range().chr().as_str(),
                entry.range().strand().as_str(),
                start,
                start + size as i32 - 1,
            )
        };
        kept.push((
            name.clone(),
            block.headers[i].clone(),
            range,
            slice.to_vec(),
        ));
    }
    if kept.is_empty() {
        return None;
    }
    let width = cend - cbeg + 1;
    let keep_col: Vec<bool> = (0..width)
        .map(|c| kept.iter().any(|(_, _, _, seq)| seq[c] != b'-'))
        .collect();
    let mut entries = Vec::new();
    let mut names = Vec::new();
    let mut headers = Vec::new();
    for (name, header, range, seq) in kept {
        let new_seq: Vec<u8> = keep_col
            .iter()
            .zip(seq.iter())
            .filter(|(k, _)| **k)
            .map(|(_, b)| *b)
            .collect();
        entries.push(FasEntry::from(&range, &new_seq));
        names.push(name);
        headers.push(header);
    }
    Some(FasBlock {
        entries,
        names,
        headers,
    })
}

/// Slice a block to the exact columns of reference positions `[beg, end]`
/// (canonical `pre_yama` overlap slicing, no flanking-gap extension).
fn slice_overlap(block: &FasBlock, ref_name: &str, beg: i32, end: i32) -> Option<FasBlock> {
    let ref_entry = find_ref_entry(block, ref_name)?;
    let cbeg = ref_pos_to_col(ref_entry, beg)?;
    let cend = ref_pos_to_col(ref_entry, end)?;
    slice_block_cols(block, ref_name, cbeg, cend, beg, end)
}

/// Slice a block to the columns covering reference positions `[beg, end]`,
/// extended over flanking insertion columns (canonical `print_part_ali_col`
/// front/tail emission).
fn slice_part(block: &FasBlock, ref_name: &str, beg: i32, end: i32) -> Option<FasBlock> {
    let ref_entry = find_ref_entry(block, ref_name)?;
    if beg > end {
        return None;
    }
    let mut cbeg = ref_pos_to_col(ref_entry, beg)?;
    let mut cend = ref_pos_to_col(ref_entry, end)?;
    while cbeg > 0 && ref_entry.seq()[cbeg - 1] == b'-' {
        cbeg -= 1;
    }
    while cend + 1 < ref_entry.seq().len() && ref_entry.seq()[cend + 1] == b'-' {
        cend += 1;
    }
    slice_block_cols(block, ref_name, cbeg, cend, beg, end)
}

/// Truncate a block to reference position `beg` and beyond (canonical
/// `keep_ali`); None when nothing remains.
fn keep_from(block: &FasBlock, ref_name: &str, beg: i32) -> Option<FasBlock> {
    let ref_entry = find_ref_entry(block, ref_name)?;
    if beg > *ref_entry.range().end() {
        return None;
    }
    slice_part(block, ref_name, beg, *ref_entry.range().end())
}

/// Columns before the reference base at `pos` (leading insertion columns of a
/// block whose reference starts at `pos`; canonical multiz gap-front output).
fn leading_insertions(block: &FasBlock, ref_name: &str, pos: i32) -> Option<FasBlock> {
    let ref_entry = find_ref_entry(block, ref_name)?;
    let col = ref_pos_to_col(ref_entry, pos)?;
    if col == 0 {
        return None;
    }
    slice_block_cols(block, ref_name, 0, col - 1, pos, pos - 1)
}

/// Columns after the reference base at `pos` (trailing insertion columns of a
/// block whose reference ends at `pos`; canonical multiz tail output).
fn trailing_insertions(block: &FasBlock, ref_name: &str, pos: i32) -> Option<FasBlock> {
    let ref_entry = find_ref_entry(block, ref_name)?;
    let col = ref_pos_to_col(ref_entry, pos)?;
    if col + 1 >= ref_entry.seq().len() {
        return None;
    }
    slice_block_cols(
        block,
        ref_name,
        col + 1,
        ref_entry.seq().len() - 1,
        pos + 1,
        pos,
    )
}

/// Slice both blocks to the shared overlap `[beg, end]` and run the block-pair
/// DP merge (canonical `pre_yama` on the overlap region).
fn merge_overlap(
    ref_name: &str,
    b1: &FasBlock,
    b2: &FasBlock,
    beg: i32,
    end: i32,
    cfg: &FasMultizConfig,
) -> anyhow::Result<Option<FasBlock>> {
    let (Some(a1), Some(a2)) = (
        slice_overlap(b1, ref_name, beg, end),
        slice_overlap(b2, ref_name, beg, end),
    ) else {
        return Ok(None);
    };
    merge_two_blocks_with_dp(ref_name, [&a1, &a2], cfg)
}

/// Canonical multiz block-stream merge (mirrors `multiz.c::multiz`): walk two
/// sorted single-coverage streams, emit non-overlapping front blocks, slice
/// out and DP-merge each overlapping region, and carry unconsumed tails into
/// the next iteration. Returns multiple output blocks in reference order.
fn merge_two_streams(
    ref_name: &str,
    list1: Vec<FasBlock>,
    list2: Vec<FasBlock>,
    cfg: &FasMultizConfig,
) -> anyhow::Result<Vec<FasBlock>> {
    let mut out: Vec<FasBlock> = Vec::new();
    let mut l1 = list1.into_iter();
    let mut l2 = list2.into_iter();
    let mut a1 = l1.next();
    let mut a2 = l2.next();

    loop {
        // Emit front blocks entirely before the other stream's current block.
        while let Some(b1) = a1.take() {
            // Blocks without a reference entry (insertion-only parts) cannot
            // overlap anything and pass straight through.
            let overlaps = match (
                ref_range(&b1, ref_name),
                a2.as_ref().and_then(|b2| ref_range(b2, ref_name)),
            ) {
                (Some((_, e1)), Some((s2, _))) => e1 >= s2,
                _ => false,
            };
            if overlaps {
                a1 = Some(b1);
                break;
            }
            if ref_size(&b1, ref_name) >= cfg.min_width as i32 {
                out.push(b1);
            }
            a1 = l1.next();
        }
        while let Some(b2) = a2.take() {
            let overlaps = match (
                ref_range(&b2, ref_name),
                a1.as_ref().and_then(|b1| ref_range(b1, ref_name)),
            ) {
                (Some((_, e2)), Some((s1, _))) => e2 >= s1,
                _ => false,
            };
            if overlaps {
                a2 = Some(b2);
                break;
            }
            if ref_size(&b2, ref_name) >= cfg.min_width as i32 {
                out.push(b2);
            }
            a2 = l2.next();
        }

        // Re-check overlap after advancing either stream: the front loops only
        // tested the previous head, so a block may now sit entirely before the
        // other stream's head (canonical multiz re-checks with `continue`).
        if let (Some(b1), Some(b2)) = (&a1, &a2) {
            let (s1, e1) = ref_range(b1, ref_name).unwrap();
            let (s2, e2) = ref_range(b2, ref_name).unwrap();
            if e1 < s2 || e2 < s1 {
                continue;
            }
        }

        match (a1.take(), a2.take()) {
            (Some(b1), Some(b2)) => {
                let (beg1, end1) = ref_range(&b1, ref_name).unwrap();
                let (beg2, end2) = ref_range(&b2, ref_name).unwrap();

                // Emit the earlier-starting block's front part before the
                // overlap.
                if beg1 < beg2 {
                    if beg2 - beg1 >= cfg.min_width as i32 {
                        if let Some(part) = slice_part(&b1, ref_name, beg1, beg2 - 1) {
                            out.push(part);
                        }
                    }
                } else if beg2 < beg1 && beg1 - beg2 >= cfg.min_width as i32 {
                    if let Some(part) = slice_part(&b2, ref_name, beg2, beg1 - 1) {
                        out.push(part);
                    }
                }

                let beg = beg1.max(beg2);
                let end = end1.min(end2);

                // Leading insertion columns before the overlap in each block
                // that starts at the overlap boundary.
                if beg == beg1 {
                    if let Some(part) = leading_insertions(&b1, ref_name, beg) {
                        out.push(part);
                    }
                }
                if beg == beg2 {
                    if let Some(part) = leading_insertions(&b2, ref_name, beg) {
                        out.push(part);
                    }
                }

                // Merge the overlap region (DP on the sliced columns).
                if let Some(merged) = merge_overlap(ref_name, &b1, &b2, beg, end, cfg)? {
                    if ref_size(&merged, ref_name) >= cfg.min_width as i32 {
                        out.push(merged);
                    }
                }

                // Carry unconsumed tails into the next iteration.
                if end1 < end2 {
                    a2 = keep_from(&b2, ref_name, end1 + 1).or_else(|| l2.next());
                } else if end2 < end1 {
                    a1 = keep_from(&b1, ref_name, end2 + 1).or_else(|| l1.next());
                }
                // Emit trailing insertion columns of the fully-consumed
                // block, then advance that stream.
                if end1 <= end2 {
                    if let Some(tail) = trailing_insertions(&b1, ref_name, end1) {
                        out.push(tail);
                    }
                    a1 = l1.next();
                }
                if end2 <= end1 {
                    if let Some(tail) = trailing_insertions(&b2, ref_name, end2) {
                        out.push(tail);
                    }
                    a2 = l2.next();
                }
            }
            (rest1, rest2) => {
                // One (or both) stream exhausted: put the survivor back and
                // rerun the front loops, which emit its remaining blocks
                // (canonical multiz `continue` on a NULL stream).
                a1 = rest1;
                a2 = rest2;
                if a1.is_none() && a2.is_none() {
                    break;
                }
            }
        }
    }

    Ok(out)
}

/// Merge all blocks of every input overlapping `window`, mirroring the
/// canonical multiz block-stream merge: each input's blocks form a sorted
/// single-coverage stream, and streams are merged pairwise in input order.
/// Returns one or more output blocks in reference order.
pub fn merge_window(
    ref_name: &str,
    window: &Window,
    blocks_per_input: &[Vec<FasBlock>],
    cfg: &FasMultizConfig,
) -> anyhow::Result<Vec<FasBlock>> {
    let mut streams: Vec<Vec<FasBlock>> = Vec::new();
    for group in blocks_per_input {
        let mut blocks: Vec<FasBlock> = group
            .iter()
            .filter(|block| match find_ref_entry(block, ref_name) {
                Some(entry) => ref_overlaps_window(entry, window),
                None => false,
            })
            .cloned()
            .collect();
        blocks.sort_by_key(|b| content_key(b, ref_name));
        streams.push(blocks);
    }

    // Progressively merge streams ordered by their first block's content key,
    // so the output does not depend on the input file order (multiz uses a
    // guide tree; pgr has none, so the content-derived order stands in).
    streams.retain(|s| !s.is_empty());
    streams.sort_by_key(|s| content_key(&s[0], ref_name));
    let mut acc: Vec<FasBlock> = Vec::new();
    for stream in streams {
        if acc.is_empty() {
            acc = stream;
            continue;
        }
        acc = merge_two_streams(ref_name, acc, stream, cfg)?;
    }
    Ok(acc)
}
