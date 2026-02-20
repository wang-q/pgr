use crate::libs::chain::sub_matrix::SubMatrix;
use crate::libs::chain::GapCalc;
use crate::libs::fmt::fas::{FasBlock, FasEntry};
use std::collections::BTreeMap;
use std::path::Path;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FasMultizMode {
    Core,
    Union,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FasMultizGapModel {
    Constant,
    Medium,
    Loose,
}

#[derive(Clone, Debug)]
pub struct FasMultizConfig {
    pub ref_name: String,
    pub radius: usize,
    pub min_width: usize,
    pub mode: FasMultizMode,
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_score: i32,
    pub gap_model: FasMultizGapModel,
    pub gap_open: Option<i32>,
    pub gap_extend: Option<i32>,
    pub score_matrix: Option<String>,
}

#[derive(Clone, Debug)]
pub struct Window {
    pub chr: String,
    pub start: u64,
    pub end: u64,
}

fn find_ref_entry<'a>(block: &'a FasBlock, ref_name: &str) -> Option<&'a FasEntry> {
    block
        .entries
        .iter()
        .zip(block.names.iter())
        .find_map(|(entry, name)| if name == ref_name { Some(entry) } else { None })
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

fn ref_overlaps_window(entry: &FasEntry, window: &Window) -> bool {
    let range = entry.range();
    if range.chr() != &window.chr {
        return false;
    }
    let start = *range.start() as u64;
    let end = *range.end() as u64;
    start < window.end && end > window.start
}

fn banded_align_refs(
    blocks: [&FasBlock; 2],
    ref_name: &str,
    cfg: &FasMultizConfig,
) -> Option<(Vec<Option<usize>>, Vec<Option<usize>>)> {
    use std::cmp::min;

    let ref_a = find_ref_entry(blocks[0], ref_name)?;
    let ref_b = find_ref_entry(blocks[1], ref_name)?;

    let sa = ref_a.seq();
    let sb = ref_b.seq();

    let n = sa.len();
    let m = sb.len();

    if n == 0 || m == 0 {
        return None;
    }

    let band = cfg.radius.max(((n as isize - m as isize).unsigned_abs()) as usize);

    let width = 2 * band + 1;
    let mut score = vec![i32::MIN; (n + 1) * width];
    let mut gap_i = vec![i32::MIN; (n + 1) * width];
    let mut gap_j = vec![i32::MIN; (n + 1) * width];
    let mut trace = vec![0i8; (n + 1) * width];

    let idx = |i: usize, j: usize| -> Option<usize> {
        let band_start = if i > band { i - band } else { 0 };
        let band_end = min(m, i + band);
        if j < band_start || j > band_end {
            None
        } else {
            let offset = j + band - i;
            Some(i * width + offset)
        }
    };

    if let Some(k) = idx(0, 0) {
        score[k] = 0;
        gap_i[k] = i32::MIN;
        gap_j[k] = i32::MIN;
        trace[k] = 0;
    } else {
        return None;
    }

    let submat = if let Some(name) = &cfg.score_matrix {
        SubMatrix::from_name(name).unwrap_or_else(|_| SubMatrix::hoxd55())
    } else {
        SubMatrix::hoxd55()
    };

    let (gap_open_pen, gap_extend_pen) = if let (Some(open), Some(extend)) = (cfg.gap_open, cfg.gap_extend) {
        let scale = cfg.match_score as f64 / 100.0;
        let open_scaled = (open as f64 * scale).round() as i32;
        let extend_scaled = (extend as f64 * scale).round() as i32;
        (-open_scaled, -extend_scaled)
    } else {
        match cfg.gap_model {
            FasMultizGapModel::Constant => (cfg.gap_score, cfg.gap_score),
            FasMultizGapModel::Medium | FasMultizGapModel::Loose => {
                let gap_calc = match cfg.gap_model {
                    FasMultizGapModel::Medium => GapCalc::medium(),
                    FasMultizGapModel::Loose => GapCalc::loose(),
                    FasMultizGapModel::Constant => unreachable!(),
                };
                let c1 = gap_calc.calc(1, 0).max(1);
                let c2 = gap_calc.calc(2, 0).max(c1 + 1);
                let open_raw = 2 * c1 - c2;
                let extend_raw = c2 - c1;
                let scale = cfg.match_score as f64 / 100.0;
                let open_scaled = (open_raw as f64 * scale).round() as i32;
                let extend_scaled = (extend_raw as f64 * scale).round() as i32;
                (-open_scaled, -extend_scaled)
            }
        }
    };

    let mut profiles: Vec<(&[u8], &[u8])> = Vec::new();
    let mut map_a: BTreeMap<&str, &FasEntry> = BTreeMap::new();
    for (entry, name) in blocks[0].entries.iter().zip(blocks[0].names.iter()) {
        map_a.insert(name.as_str(), entry);
    }
    for (entry, name) in blocks[1].entries.iter().zip(blocks[1].names.iter()) {
        if let Some(ea) = map_a.get(name.as_str()) {
            profiles.push((ea.seq(), entry.seq()));
        }
    }

    for i in 0..=n {
        let band_start = if i > band { i - band } else { 0 };
        let band_end = min(m, i + band);
        for j in band_start..=band_end {
            let k = match idx(i, j) {
                Some(v) => v,
                None => continue,
            };
            if i == 0 && j == 0 {
                continue;
            }

            let mut best = i32::MIN;
            let mut bt = 0i8;

            let mut m_val = i32::MIN;
            if i > 0 && j > 0 {
                if let Some(pk) = idx(i - 1, j - 1) {
                    let mut s = 0;
                    for (pa, pb) in &profiles {
                        let ba = pa[i - 1];
                        let bb = pb[j - 1];
                        if ba == b'-' && bb == b'-' {
                            continue;
                        } else if ba == b'-' || bb == b'-' {
                            s += gap_open_pen + gap_extend_pen;
                        } else {
                            let raw = submat.get_score(ba as char, bb as char);
                            s += raw / 50;
                        }
                    }
                    m_val = score[pk].saturating_add(s);
                    best = m_val;
                    bt = 1;
                }
            }

            let mut gi_val = i32::MIN;
            if i > 0 {
                if let Some(pk_score) = idx(i - 1, j) {
                    let from_match = score[pk_score]
                        .saturating_add(gap_open_pen)
                        .saturating_add(gap_extend_pen);
                    let from_gap = gap_i[pk_score].saturating_add(gap_extend_pen);
                    gi_val = from_match.max(from_gap);
                    if gi_val > best {
                        best = gi_val;
                        bt = 2;
                    }
                }
            }

            let mut gj_val = i32::MIN;
            if j > 0 {
                if let Some(pk_score) = idx(i, j - 1) {
                    let from_match = score[pk_score]
                        .saturating_add(gap_open_pen)
                        .saturating_add(gap_extend_pen);
                    let from_gap = gap_j[pk_score].saturating_add(gap_extend_pen);
                    gj_val = from_match.max(from_gap);
                    if gj_val > best {
                        best = gj_val;
                        bt = 3;
                    }
                }
            }

            score[k] = best;
            gap_i[k] = gi_val;
            gap_j[k] = gj_val;
            trace[k] = bt;
        }
    }

    let mut i = n;
    let mut j = m;

    if idx(i, j).is_none() {
        return None;
    }

    let mut map_a = Vec::new();
    let mut map_b = Vec::new();

    while i > 0 || j > 0 {
        let k = match idx(i, j) {
            Some(v) => v,
            None => break,
        };
        let bt = trace[k];
        if bt == 1 {
            if i == 0 || j == 0 {
                break;
            }
            let pi = i - 1;
            let pj = j - 1;
            map_a.push(Some(pi));
            map_b.push(Some(pj));
            i -= 1;
            j -= 1;
        } else if bt == 2 {
            if i == 0 {
                break;
            }
            let pi = i - 1;
            map_a.push(Some(pi));
            map_b.push(None);
            i -= 1;
        } else if bt == 3 {
            if j == 0 {
                break;
            }
            let pj = j - 1;
            map_a.push(None);
            map_b.push(Some(pj));
            j -= 1;
        } else {
            break;
        }
    }

    map_a.reverse();
    map_b.reverse();

    if map_a.len() != map_b.len() || map_a.is_empty() {
        return None;
    }

    Some((map_a, map_b))
}

fn merge_two_blocks_with_dp(
    ref_name: &str,
    blocks: [&FasBlock; 2],
    cfg: &FasMultizConfig,
) -> Option<FasBlock> {
    let ref_a = find_ref_entry(blocks[0], ref_name)?;
    let ref_b = find_ref_entry(blocks[1], ref_name)?;

    if !ungapped_equal(ref_a, ref_b) {
        return None;
    }

    let (map_a, map_b) = banded_align_refs(blocks, ref_name, cfg)?;

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
            let chosen = if group[0].is_some() { group[0] } else { group[1] }.unwrap();
            chosen.range().clone()
        };

        let entry = FasEntry::from(&range, &seq);
        let header = format!("{}", range);

        entries.push(entry);
        names.push(name.clone());
        headers.push(header);
    }

    if entries.is_empty() {
        None
    } else {
        Some(FasBlock {
            entries,
            names,
            headers,
        })
    }
}

fn merge_blocks_with_dp(
    ref_name: &str,
    blocks: &[&FasBlock],
    cfg: &FasMultizConfig,
) -> Option<FasBlock> {
    if blocks.len() < 2 {
        return None;
    }

    let mut acc = merge_two_blocks_with_dp(ref_name, [blocks[0], blocks[1]], cfg)?;

    if blocks.len() == 2 {
        return Some(acc);
    }

    match cfg.mode {
        FasMultizMode::Core => {
            for &block in &blocks[2..] {
                acc = merge_two_blocks_with_dp(ref_name, [&acc, block], cfg)?;
            }
        }
        FasMultizMode::Union => {
            for &block in &blocks[2..] {
                if let Some(next) = merge_two_blocks_with_dp(ref_name, [&acc, block], cfg) {
                    acc = next;
                }
            }
        }
    }

    Some(acc)
}

pub fn merge_window(
    ref_name: &str,
    window: &Window,
    blocks_per_input: &[Vec<FasBlock>],
    cfg: &FasMultizConfig,
) -> Option<FasBlock> {
    if blocks_per_input.is_empty() {
        return None;
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
                    return None;
                }
            }
        }
    }

    if blocks.is_empty() {
        return None;
    }

    if blocks.len() >= 2 {
        if let Some(block) = merge_blocks_with_dp(ref_name, &blocks, cfg) {
            return Some(block);
        }
    }

    let template = blocks[0];
    let ref_entry = find_ref_entry(template, ref_name)?;

    for block in &blocks[1..] {
        let other_ref = find_ref_entry(block, ref_name)?;
        if !entry_seq_equal(ref_entry, other_ref) {
            return None;
        }
    }

    let ref_range = ref_entry.range().clone();

    let n = blocks.len();
    let mut species_map: BTreeMap<String, Vec<Option<&FasEntry>>> = BTreeMap::new();

    for (i, block) in blocks.iter().enumerate() {
        for (entry, name) in block.entries.iter().zip(block.names.iter()) {
            let v = species_map.entry(name.clone()).or_insert_with(|| vec![None; n]);
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

        let seq = chosen.seq().clone();
        let entry = FasEntry::from(&range, &seq);
        let header = format!("{}", range);

        entries.push(entry);
        names.push(name.clone());
        headers.push(header);
    }

    if entries.is_empty() {
        None
    } else {
        Some(FasBlock {
            entries,
            names,
            headers,
        })
    }
}

pub fn merge_fas_files(
    ref_name: &str,
    infiles: &[impl AsRef<Path>],
    windows: &[Window],
    cfg: &FasMultizConfig,
) -> anyhow::Result<Vec<FasBlock>> {
    use std::fs::File;
    use std::io::BufReader;

    let mut blocks_per_input: Vec<Vec<FasBlock>> = Vec::new();

    for infile in infiles {
        let file = File::open(infile)?;
        let mut reader = BufReader::new(file);
        let mut blocks = Vec::new();

        loop {
            match crate::libs::fmt::fas::next_fas_block(&mut reader) {
                Ok(block) => blocks.push(block),
                Err(e) => {
                    if e.to_string() == "EOF" {
                        break;
                    } else {
                        return Err(e.into());
                    }
                }
            }
        }

        blocks_per_input.push(blocks);
    }

    if windows.is_empty() {
        return Ok(Vec::new());
    }

    let mut merged_blocks = Vec::new();
    for window in windows {
        if let Some(block) = merge_window(ref_name, window, &blocks_per_input, cfg) {
            merged_blocks.push(block);
        }
    }

    Ok(merged_blocks)
}

fn derive_windows_from_blocks(
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

pub fn merge_fas_files_auto_windows(
    ref_name: &str,
    infiles: &[impl AsRef<Path>],
    cfg: &FasMultizConfig,
) -> anyhow::Result<Vec<FasBlock>> {
    use std::fs::File;
    use std::io::BufReader;

    let mut blocks_per_input: Vec<Vec<FasBlock>> = Vec::new();

    for infile in infiles {
        let file = File::open(infile)?;
        let mut reader = BufReader::new(file);
        let mut blocks = Vec::new();

        loop {
            match crate::libs::fmt::fas::next_fas_block(&mut reader) {
                Ok(block) => blocks.push(block),
                Err(e) => {
                    if e.to_string() == "EOF" {
                        break;
                    } else {
                        return Err(e.into());
                    }
                }
            }
        }

        blocks_per_input.push(blocks);
    }

    let windows = derive_windows_from_blocks(ref_name, &blocks_per_input, cfg);
    if windows.is_empty() {
        return Ok(Vec::new());
    }

    let mut merged_blocks = Vec::new();
    for window in &windows {
        if let Some(block) = merge_window(ref_name, window, &blocks_per_input, cfg) {
            merged_blocks.push(block);
        }
    }

    Ok(merged_blocks)
}

#[cfg(test)]
mod tests {
    use super::*;
    use intspan::Range;

    fn make_entry(name: &str, start: i32, end: i32, seq: &str) -> (FasEntry, String, String) {
        let range = Range::from(name, start, end);
        let entry = FasEntry::from(&range, seq.as_bytes());
        let header = format!("{}", range);
        (entry, name.to_string(), header)
    }

    fn make_block(entries: Vec<(FasEntry, String, String)>) -> FasBlock {
        let mut es = Vec::new();
        let mut names = Vec::new();
        let mut headers = Vec::new();
        for (e, n, h) in entries {
            es.push(e);
            names.push(n);
            headers.push(h);
        }
        FasBlock {
            entries: es,
            names,
            headers,
        }
    }

    fn default_config(mode: FasMultizMode) -> FasMultizConfig {
        FasMultizConfig {
            ref_name: "ref".to_string(),
            radius: 5,
            min_width: 1,
            mode,
            match_score: 2,
            mismatch_score: -1,
            gap_score: -2,
            gap_model: FasMultizGapModel::Medium,
            gap_open: None,
            gap_extend: None,
            score_matrix: None,
        }
    }

    #[test]
    fn merge_window_core_requires_blocks_in_all_inputs() {
        let (ref_entry, ref_name, ref_header) = make_entry("ref", 1, 4, "ACGT");
        let (a_entry, a_name, a_header) = make_entry("A", 1, 4, "ACGT");
        let block1 = make_block(vec![
            (ref_entry.clone(), ref_name.clone(), ref_header.clone()),
            (a_entry, a_name, a_header),
        ]);

        let blocks_per_input = vec![vec![block1], Vec::new()];

        let cfg = default_config(FasMultizMode::Union);
        let window = Window {
            chr: "chr1".to_string(),
            start: 1,
            end: 4,
        };

        let merged = merge_window("ref", &window, &blocks_per_input, &cfg);
        assert!(merged.is_none());
    }

    #[test]
    fn merge_window_core_requires_overlap_with_window() {
        let (ref_entry, ref_name, ref_header) = make_entry("ref", 1, 4, "ACGT");
        let (a_entry, a_name, a_header) = make_entry("A", 1, 4, "ACGT");
        let block1 = make_block(vec![
            (ref_entry.clone(), ref_name.clone(), ref_header.clone()),
            (a_entry, a_name, a_header),
        ]);

        let blocks_per_input = vec![vec![block1]];

        let cfg = default_config(FasMultizMode::Core);
        let window = Window {
            chr: ref_entry.range().chr().to_string(),
            start: 10,
            end: 20,
        };

        let merged = merge_window("ref", &window, &blocks_per_input, &cfg);
        assert!(merged.is_none());
    }

    #[test]
    fn merge_window_union_allows_missing_blocks() {
        let (ref_entry, ref_name, ref_header) = make_entry("ref", 1, 4, "ACGT");
        let (a_entry, a_name, a_header) = make_entry("A", 1, 4, "ACGT");
        let block1 = make_block(vec![
            (ref_entry.clone(), ref_name.clone(), ref_header.clone()),
            (a_entry, a_name, a_header),
        ]);

        let blocks_per_input = vec![vec![block1], Vec::new()];

        let cfg = default_config(FasMultizMode::Union);
        let window = Window {
            chr: ref_entry.range().chr().to_string(),
            start: *ref_entry.range().start() as u64,
            end: *ref_entry.range().end() as u64,
        };

        let merged = merge_window("ref", &window, &blocks_per_input, &cfg).unwrap();
        assert_eq!(merged.names.len(), 2);
        assert_eq!(merged.names[0], "ref");
        assert_eq!(merged.names[1], "A");
    }

    #[test]
    fn merge_window_core_species_intersection() {
        let (ref_entry1, ref_name1, ref_header1) = make_entry("ref", 1, 4, "ACGT");
        let (a_entry1, a_name1, a_header1) = make_entry("A", 1, 4, "ACGT");
        let (b_entry1, b_name1, b_header1) = make_entry("B", 1, 4, "ACGT");
        let block1 = make_block(vec![
            (ref_entry1.clone(), ref_name1.clone(), ref_header1.clone()),
            (a_entry1, a_name1, a_header1),
            (b_entry1, b_name1, b_header1),
        ]);

        let (ref_entry2, ref_name2, ref_header2) = make_entry("ref", 1, 4, "ACGT");
        let (a_entry2, a_name2, a_header2) = make_entry("A", 1, 4, "ACGT");
        let (c_entry2, c_name2, c_header2) = make_entry("C", 1, 4, "ACGT");
        let block2 = make_block(vec![
            (ref_entry2.clone(), ref_name2.clone(), ref_header2.clone()),
            (a_entry2, a_name2, a_header2),
            (c_entry2, c_name2, c_header2),
        ]);

        let blocks_per_input = vec![vec![block1], vec![block2]];

        let cfg = default_config(FasMultizMode::Core);
        let window = Window {
            chr: ref_entry1.range().chr().to_string(),
            start: *ref_entry1.range().start() as u64,
            end: *ref_entry1.range().end() as u64,
        };

        let merged = merge_window("ref", &window, &blocks_per_input, &cfg).unwrap();

        let names: Vec<String> = merged.names.clone();
        assert!(names.contains(&"ref".to_string()));
        assert!(names.contains(&"A".to_string()));
        assert_eq!(names.len(), 2);
    }

    #[test]
    fn merge_window_union_species_union() {
        let (ref_entry1, ref_name1, ref_header1) = make_entry("ref", 1, 4, "ACGT");
        let (a_entry1, a_name1, a_header1) = make_entry("A", 1, 4, "ACGT");
        let (b_entry1, b_name1, b_header1) = make_entry("B", 1, 4, "ACGT");
        let block1 = make_block(vec![
            (ref_entry1.clone(), ref_name1.clone(), ref_header1.clone()),
            (a_entry1, a_name1, a_header1),
            (b_entry1, b_name1, b_header1),
        ]);

        let (ref_entry2, ref_name2, ref_header2) = make_entry("ref", 1, 4, "ACGT");
        let (a_entry2, a_name2, a_header2) = make_entry("A", 1, 4, "ACGT");
        let (c_entry2, c_name2, c_header2) = make_entry("C", 1, 4, "ACGT");
        let block2 = make_block(vec![
            (ref_entry2.clone(), ref_name2.clone(), ref_header2.clone()),
            (a_entry2, a_name2, a_header2),
            (c_entry2, c_name2, c_header2),
        ]);

        let blocks_per_input = vec![vec![block1], vec![block2]];

        let cfg = default_config(FasMultizMode::Union);
        let window = Window {
            chr: ref_entry1.range().chr().to_string(),
            start: *ref_entry1.range().start() as u64,
            end: *ref_entry1.range().end() as u64,
        };

        let merged = merge_window("ref", &window, &blocks_per_input, &cfg).unwrap();

        let mut names: Vec<String> = merged.names.clone();
        names.sort();
        assert_eq!(names, vec!["A".to_string(), "B".to_string(), "C".to_string(), "ref".to_string()]);
    }

    #[test]
    fn merge_window_mismatched_reference_returns_none() {
        let (ref_entry1, ref_name1, ref_header1) = make_entry("ref", 1, 4, "ACGT");
        let (a_entry1, a_name1, a_header1) = make_entry("A", 1, 4, "ACGT");
        let block1 = make_block(vec![
            (ref_entry1, ref_name1, ref_header1),
            (a_entry1, a_name1, a_header1),
        ]);

        let (ref_entry2, ref_name2, ref_header2) = make_entry("ref", 1, 4, "AGGT");
        let (a_entry2, a_name2, a_header2) = make_entry("A", 1, 4, "AGGT");
        let block2 = make_block(vec![
            (ref_entry2, ref_name2, ref_header2),
            (a_entry2, a_name2, a_header2),
        ]);

        let blocks_per_input = vec![vec![block1], vec![block2]];

        let cfg = default_config(FasMultizMode::Union);
        let window = Window {
            chr: "ref".to_string(),
            start: 1,
            end: 4,
        };

        let merged = merge_window("ref", &window, &blocks_per_input, &cfg);
        assert!(merged.is_none());
    }

    #[test]
    fn merge_fas_files_multiple_windows() {
        use std::fs::File;
        use std::io::Write;
        use intspan::Range;

        let dir = std::env::temp_dir();
        let path1 = dir.join("pgr_fas_multiz_test1.fas");
        let path2 = dir.join("pgr_fas_multiz_test2.fas");

        {
            let mut f1 = File::create(&path1).unwrap();
            writeln!(f1, ">ref(+):1-4|species=ref").unwrap();
            writeln!(f1, "ACGT").unwrap();
            writeln!(f1, ">A(+):1-4|species=A").unwrap();
            writeln!(f1, "ACGT").unwrap();
            writeln!(f1).unwrap();
            writeln!(f1, ">ref(+):21-24|species=ref").unwrap();
            writeln!(f1, "ACGT").unwrap();
            writeln!(f1, ">A(+):21-24|species=A").unwrap();
            writeln!(f1, "ACGT").unwrap();
            writeln!(f1).unwrap();

            let mut f2 = File::create(&path2).unwrap();
            writeln!(f2, ">ref(+):1-4|species=ref").unwrap();
            writeln!(f2, "ACGT").unwrap();
            writeln!(f2, ">B(+):1-4|species=B").unwrap();
            writeln!(f2, "ACGT").unwrap();
            writeln!(f2).unwrap();
            writeln!(f2, ">ref(+):21-24|species=ref").unwrap();
            writeln!(f2, "ACGT").unwrap();
            writeln!(f2, ">B(+):21-24|species=B").unwrap();
            writeln!(f2, "ACGT").unwrap();
            writeln!(f2).unwrap();
        }

        let r1 = Range::from_str("ref(+):1-4|species=ref");
        let r2 = Range::from_str("ref(+):21-24|species=ref");

        let windows = vec![
            Window {
                chr: r1.chr().to_string(),
                start: *r1.start() as u64,
                end: *r1.end() as u64,
            },
            Window {
                chr: r2.chr().to_string(),
                start: *r2.start() as u64,
                end: *r2.end() as u64,
            },
        ];

        let ref_name = "ref".to_string();
        let mut cfg = default_config(FasMultizMode::Union);
        cfg.ref_name = ref_name.clone();

        let merged =
            merge_fas_files(&ref_name, &[&path1, &path2], &windows, &cfg).expect("merge_fas_files");

        assert_eq!(merged.len(), 2);
        for block in merged {
            let mut names = block.names.clone();
            names.sort();
            assert_eq!(names, vec!["A".to_string(), "B".to_string(), "ref".to_string()]);
        }
    }

    #[test]
    fn merge_fas_files_auto_windows_matches_explicit() {
        use std::fs::File;
        use std::io::Write;
        use intspan::Range;

        let dir = std::env::temp_dir();
        let path1 = dir.join("pgr_fas_multiz_auto_test1.fas");
        let path2 = dir.join("pgr_fas_multiz_auto_test2.fas");

        {
            let mut f1 = File::create(&path1).unwrap();
            writeln!(f1, ">ref(+):1-4|species=ref").unwrap();
            writeln!(f1, "ACGT").unwrap();
            writeln!(f1, ">A(+):1-4|species=A").unwrap();
            writeln!(f1, "ACGT").unwrap();
            writeln!(f1).unwrap();
            writeln!(f1, ">ref(+):21-24|species=ref").unwrap();
            writeln!(f1, "ACGT").unwrap();
            writeln!(f1, ">A(+):21-24|species=A").unwrap();
            writeln!(f1, "ACGT").unwrap();
            writeln!(f1).unwrap();

            let mut f2 = File::create(&path2).unwrap();
            writeln!(f2, ">ref(+):1-4|species=ref").unwrap();
            writeln!(f2, "ACGT").unwrap();
            writeln!(f2, ">B(+):1-4|species=B").unwrap();
            writeln!(f2, "ACGT").unwrap();
            writeln!(f2).unwrap();
            writeln!(f2, ">ref(+):21-24|species=ref").unwrap();
            writeln!(f2, "ACGT").unwrap();
            writeln!(f2, ">B(+):21-24|species=B").unwrap();
            writeln!(f2, "ACGT").unwrap();
            writeln!(f2).unwrap();
        }

        let r1 = Range::from_str("ref(+):1-4|species=ref");
        let r2 = Range::from_str("ref(+):21-24|species=ref");

        let windows = vec![
            Window {
                chr: r1.chr().to_string(),
                start: *r1.start() as u64,
                end: *r1.end() as u64,
            },
            Window {
                chr: r2.chr().to_string(),
                start: *r2.start() as u64,
                end: *r2.end() as u64,
            },
        ];

        let ref_name = "ref".to_string();
        let mut cfg = default_config(FasMultizMode::Union);
        cfg.ref_name = ref_name.clone();

        let merged_explicit =
            merge_fas_files(&ref_name, &[&path1, &path2], &windows, &cfg).expect("merge_fas_files");
        let merged_auto = merge_fas_files_auto_windows(&ref_name, &[&path1, &path2], &cfg)
            .expect("merge_fas_files_auto_windows");

        assert_eq!(merged_explicit.len(), merged_auto.len());
        for (a, b) in merged_explicit.iter().zip(merged_auto.iter()) {
            assert_eq!(a.names, b.names);
        }
    }

    #[test]
    fn merge_window_multi_input_dp_progressive() {
        let (ref_entry1, ref_name1, ref_header1) = make_entry("ref", 1, 4, "AC--GT");
        let (a_entry1, a_name1, a_header1) = make_entry("A", 1, 4, "AC--GT");
        let block1 = make_block(vec![
            (ref_entry1.clone(), ref_name1.clone(), ref_header1.clone()),
            (a_entry1, a_name1, a_header1),
        ]);

        let (ref_entry2, ref_name2, ref_header2) = make_entry("ref", 1, 4, "A-C-GT");
        let (b_entry2, b_name2, b_header2) = make_entry("B", 1, 4, "A-C-GT");
        let block2 = make_block(vec![
            (ref_entry2.clone(), ref_name2.clone(), ref_header2.clone()),
            (b_entry2, b_name2, b_header2),
        ]);

        let (ref_entry3, ref_name3, ref_header3) = make_entry("ref", 1, 4, "ACG-T-");
        let (c_entry3, c_name3, c_header3) = make_entry("C", 1, 4, "ACG-T-");
        let block3 = make_block(vec![
            (ref_entry3, ref_name3, ref_header3),
            (c_entry3, c_name3, c_header3),
        ]);

        let blocks_per_input = vec![vec![block1], vec![block2], vec![block3]];

        let cfg = default_config(FasMultizMode::Union);
        let window = Window {
            chr: "ref".to_string(),
            start: 1,
            end: 4,
        };

        let merged = merge_window("ref", &window, &blocks_per_input, &cfg).unwrap();
        let mut names = merged.names.clone();
        names.sort();
        assert_eq!(
            names,
            vec!["A".to_string(), "B".to_string(), "C".to_string(), "ref".to_string()]
        );
    }
}
