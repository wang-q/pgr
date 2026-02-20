use crate::libs::fmt::fas::{FasBlock, FasEntry};
use std::collections::BTreeMap;
use std::path::Path;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FasMultizMode {
    Core,
    Union,
}

#[derive(Clone, Debug)]
pub struct FasMultizConfig {
    pub ref_name: String,
    pub radius: usize,
    pub min_width: usize,
    pub mode: FasMultizMode,
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

fn ref_overlaps_window(entry: &FasEntry, window: &Window) -> bool {
    let range = entry.range();
    if range.chr() != &window.chr {
        return false;
    }
    let start = *range.start() as u64;
    let end = *range.end() as u64;
    start < window.end && end > window.start
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
}
