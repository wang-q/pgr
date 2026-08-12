use super::*;
use crate::libs::ds::Range;
use crate::libs::fas_multiz::windows::derive_windows_from_blocks;

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

fn default_config() -> FasMultizConfig {
    FasMultizConfig {
        ref_name: "ref".to_string(),
        radius: 5,
        min_width: 1,
    }
}

#[test]
fn derive_windows_inverted_reference_range_no_panic() {
    // A malformed but parseable header like `>ref.chr1(+):100-1` yields an
    // inverted reference range (start > end). Deriving windows from it used to
    // underflow `width = e - s` (debug panic, release wrap). It must be
    // skipped without panicking.
    let (ref_entry, ref_name, ref_header) = make_entry("ref", 100, 1, "ACGT");
    let block = make_block(vec![(ref_entry, ref_name, ref_header)]);
    let blocks_per_input = vec![vec![block]];

    let cfg = default_config();
    let windows = derive_windows_from_blocks("ref", &blocks_per_input, &cfg);
    assert!(windows.is_empty());
}

#[test]
fn ref_overlaps_window_includes_boundary_positions() {
    // Both the block range and the window are 1-based and inclusive, so a
    // single-base block sitting exactly on the window's start/end is inside.
    let (first, _, _) = make_entry("ref", 20, 20, "A");
    let window = Window {
        chr: first.range().chr().to_string(),
        start: 10,
        end: 20,
    };
    // base 20 is the window's inclusive end -> must overlap.
    assert!(ref_overlaps_window(&first, &window));
    // base 10 is the window's inclusive start -> must overlap.
    let (entry_at_start, _, _) = make_entry("ref", 10, 10, "A");
    assert!(ref_overlaps_window(&entry_at_start, &window));
    // just outside -> must not overlap.
    let (entry_before, _, _) = make_entry("ref", 9, 9, "A");
    assert!(!ref_overlaps_window(&entry_before, &window));
    let (entry_after, _, _) = make_entry("ref", 21, 21, "A");
    assert!(!ref_overlaps_window(&entry_after, &window));
}

#[test]
fn merge_window_without_coverage_returns_none() {
    let (ref_entry, ref_name, ref_header) = make_entry("ref", 1, 4, "ACGT");
    let (a_entry, a_name, a_header) = make_entry("A", 1, 4, "ACGT");
    let block1 = make_block(vec![
        (ref_entry.clone(), ref_name.clone(), ref_header.clone()),
        (a_entry, a_name, a_header),
    ]);

    let blocks_per_input = vec![vec![block1]];

    let cfg = default_config();
    let window = Window {
        chr: ref_entry.range().chr().to_string(),
        start: 10,
        end: 20,
    };

    let merged = merge_window("ref", &window, &blocks_per_input, &cfg).unwrap();
    assert!(merged.is_empty());
}

#[test]
fn merge_window_skips_missing_inputs() {
    let (ref_entry, ref_name, ref_header) = make_entry("ref", 1, 4, "ACGT");
    let (a_entry, a_name, a_header) = make_entry("A", 1, 4, "ACGT");
    let block1 = make_block(vec![
        (ref_entry.clone(), ref_name.clone(), ref_header.clone()),
        (a_entry, a_name, a_header),
    ]);

    let blocks_per_input = vec![vec![block1], Vec::new()];

    let cfg = default_config();
    let window = Window {
        chr: ref_entry.range().chr().to_string(),
        start: *ref_entry.range().start() as u64,
        end: *ref_entry.range().end() as u64,
    };

    let merged = merge_window("ref", &window, &blocks_per_input, &cfg).unwrap();
    let merged = merged.into_iter().next().expect("one merged block");
    assert_eq!(merged.names.len(), 2);
    assert_eq!(merged.names[0], "ref");
    assert_eq!(merged.names[1], "A");
}

#[test]
fn merge_window_keeps_species_union() {
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

    let cfg = default_config();
    let window = Window {
        chr: ref_entry1.range().chr().to_string(),
        start: *ref_entry1.range().start() as u64,
        end: *ref_entry1.range().end() as u64,
    };

    let merged = merge_window("ref", &window, &blocks_per_input, &cfg).unwrap();
    let merged = merged.into_iter().next().expect("one merged block");

    let mut names: Vec<String> = merged.names.clone();
    names.sort();
    assert_eq!(
        names,
        vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "ref".to_string()
        ]
    );
}

#[test]
fn merge_window_mismatched_reference_splices_at_crossover() {
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

    let cfg = default_config();
    let window = Window {
        chr: "ref".to_string(),
        start: 1,
        end: 4,
    };

    let merged = merge_window("ref", &window, &blocks_per_input, &cfg).unwrap();
    // Mismatched reference sequences are spliced at the best crossover point
    // instead of being dropped; the merged block keeps both species.
    let block = merged
        .into_iter()
        .next()
        .expect("mismatched refs should merge via crossover");
    assert!(block.names.contains(&"ref".to_string()));
    assert!(block.names.contains(&"A".to_string()));
}

#[test]
fn merge_window_mismatched_reference_no_shared_species_returns_none() {
    let (ref_entry1, ref_name1, ref_header1) = make_entry("ref", 1, 4, "ACGT");
    let (a_entry1, a_name1, a_header1) = make_entry("A", 1, 4, "ACGT");
    let block1 = make_block(vec![
        (ref_entry1, ref_name1, ref_header1),
        (a_entry1, a_name1, a_header1),
    ]);

    let (ref_entry2, ref_name2, ref_header2) = make_entry("ref", 1, 4, "AGGT");
    let (b_entry2, b_name2, b_header2) = make_entry("B", 1, 4, "AGGT");
    let block2 = make_block(vec![
        (ref_entry2, ref_name2, ref_header2),
        (b_entry2, b_name2, b_header2),
    ]);

    let blocks_per_input = vec![vec![block1], vec![block2]];

    let cfg = default_config();
    let window = Window {
        chr: "ref".to_string(),
        start: 1,
        end: 4,
    };

    let merged = merge_window("ref", &window, &blocks_per_input, &cfg).unwrap();
    // Without a shared non-reference species there is nothing to score the
    // crossover with, so the merge is still refused.
    assert!(merged.is_empty());
}

#[test]
fn merge_window_conflicting_refs_keeps_single_block_species() {
    // When two blocks conflict in their reference sequence and are spliced at
    // a crossover, a species present in only one block must still be carried
    // across the whole merged output. Regression for a bug where the half of
    // such a species' sequence on the other side of the crossover collapsed to
    // gaps (silent data loss).
    let (ref_entry1, ref_name1, ref_header1) = make_entry("ref", 1, 4, "ACGT");
    let (a_entry1, a_name1, a_header1) = make_entry("A", 1, 4, "ACGT");
    let (b_entry1, b_name1, b_header1) = make_entry("B", 1, 4, "ACGT");
    let block1 = make_block(vec![
        (ref_entry1.clone(), ref_name1.clone(), ref_header1.clone()),
        (a_entry1, a_name1, a_header1),
        (b_entry1, b_name1, b_header1),
    ]);

    let (ref_entry2, ref_name2, ref_header2) = make_entry("ref", 1, 4, "AGGT");
    let (a_entry2, a_name2, a_header2) = make_entry("A", 1, 4, "AGGT");
    let (c_entry2, c_name2, c_header2) = make_entry("C", 1, 4, "AGGT");
    let block2 = make_block(vec![
        (ref_entry2.clone(), ref_name2.clone(), ref_header2.clone()),
        (a_entry2, a_name2, a_header2),
        (c_entry2, c_name2, c_header2),
    ]);

    let blocks_per_input = vec![vec![block1], vec![block2]];

    let cfg = default_config();
    let window = Window {
        chr: "ref".to_string(),
        start: 1,
        end: 4,
    };

    let merged = merge_window("ref", &window, &blocks_per_input, &cfg).unwrap();
    let block = merged
        .into_iter()
        .next()
        .expect("conflicting refs should merge via crossover");
    assert!(block.names.contains(&"B".to_string()));
    assert!(block.names.contains(&"C".to_string()));

    // Species B (block A only) and C (block B only) must keep their full
    // sequences; neither may collapse to all gaps on one side of the splice.
    let idx_b = block.names.iter().position(|n| n == "B").unwrap();
    let idx_c = block.names.iter().position(|n| n == "C").unwrap();
    let seq_b = block.entries[idx_b].seq();
    let seq_c = block.entries[idx_c].seq();
    assert_eq!(seq_b.len(), seq_c.len());
    assert!(
        seq_b.iter().any(|b| *b != b'-'),
        "B truncated to gaps: {:?}",
        seq_b
    );
    assert!(
        seq_c.iter().any(|b| *b != b'-'),
        "C truncated to gaps: {:?}",
        seq_c
    );
}

#[test]
fn merge_fas_files_multiple_windows() {
    use crate::libs::ds::Range;
    use std::fs::File;
    use std::io::Write;

    let dir = std::env::temp_dir();
    let path1 = dir.join("pgr_fas_multiz_test1.fas");
    let path2 = dir.join("pgr_fas_multiz_test2.fas");

    {
        let mut f1 = File::create(&path1).unwrap();
        writeln!(f1, ">ref.chr1(+):1-4").unwrap();
        writeln!(f1, "ACGT").unwrap();
        writeln!(f1, ">A.chr1(+):1-4").unwrap();
        writeln!(f1, "ACGT").unwrap();
        writeln!(f1).unwrap();
        writeln!(f1, ">ref.chr1(+):21-24").unwrap();
        writeln!(f1, "ACGT").unwrap();
        writeln!(f1, ">A.chr1(+):21-24").unwrap();
        writeln!(f1, "ACGT").unwrap();
        writeln!(f1).unwrap();

        let mut f2 = File::create(&path2).unwrap();
        writeln!(f2, ">ref.chr1(+):1-4").unwrap();
        writeln!(f2, "ACGT").unwrap();
        writeln!(f2, ">B.chr1(+):1-4").unwrap();
        writeln!(f2, "ACGT").unwrap();
        writeln!(f2).unwrap();
        writeln!(f2, ">ref.chr1(+):21-24").unwrap();
        writeln!(f2, "ACGT").unwrap();
        writeln!(f2, ">B.chr1(+):21-24").unwrap();
        writeln!(f2, "ACGT").unwrap();
        writeln!(f2).unwrap();
    }

    let r1 = Range::from_str("ref.chr1(+):1-4");
    let r2 = Range::from_str("ref.chr1(+):21-24");

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
    let mut cfg = default_config();
    cfg.ref_name = ref_name.clone();

    let merged =
        merge_fas_files(&ref_name, &[&path1, &path2], &windows, &cfg).expect("merge_fas_files");

    assert_eq!(merged.len(), 2);
    for block in merged {
        let mut names = block.names.clone();
        names.sort();
        assert_eq!(
            names,
            vec!["A".to_string(), "B".to_string(), "ref".to_string()]
        );
    }
}

#[test]
fn merge_fas_files_auto_windows_matches_explicit() {
    use crate::libs::ds::Range;
    use std::fs::File;
    use std::io::Write;

    let dir = std::env::temp_dir();
    let path1 = dir.join("pgr_fas_multiz_auto_test1.fas");
    let path2 = dir.join("pgr_fas_multiz_auto_test2.fas");

    {
        let mut f1 = File::create(&path1).unwrap();
        writeln!(f1, ">ref.chr1(+):1-4").unwrap();
        writeln!(f1, "ACGT").unwrap();
        writeln!(f1, ">A.chr1(+):1-4").unwrap();
        writeln!(f1, "ACGT").unwrap();
        writeln!(f1).unwrap();
        writeln!(f1, ">ref.chr1(+):21-24").unwrap();
        writeln!(f1, "ACGT").unwrap();
        writeln!(f1, ">A.chr1(+):21-24").unwrap();
        writeln!(f1, "ACGT").unwrap();
        writeln!(f1).unwrap();

        let mut f2 = File::create(&path2).unwrap();
        writeln!(f2, ">ref.chr1(+):1-4").unwrap();
        writeln!(f2, "ACGT").unwrap();
        writeln!(f2, ">B.chr1(+):1-4").unwrap();
        writeln!(f2, "ACGT").unwrap();
        writeln!(f2).unwrap();
        writeln!(f2, ">ref.chr1(+):21-24").unwrap();
        writeln!(f2, "ACGT").unwrap();
        writeln!(f2, ">B.chr1(+):21-24").unwrap();
        writeln!(f2, "ACGT").unwrap();
        writeln!(f2).unwrap();
    }

    let r1 = Range::from_str("ref.chr1(+):1-4");
    let r2 = Range::from_str("ref.chr1(+):21-24");

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
    let mut cfg = default_config();
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

    let cfg = default_config();
    let window = Window {
        chr: "ref".to_string(),
        start: 1,
        end: 4,
    };

    let merged = merge_window("ref", &window, &blocks_per_input, &cfg)
        .unwrap()
        .into_iter()
        .next()
        .expect("one merged block");
    let mut names = merged.names.clone();
    names.sort();
    assert_eq!(
        names,
        vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "ref".to_string()
        ]
    );
}

#[test]
fn merge_window_output_independent_of_input_order() {
    // Same three blocks fed in every input-file permutation must produce an
    // identical merged block (names and sequences). The progressive merge
    // order is derived from block contents, never from input order.
    let perms = [
        [0usize, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];

    let build_blocks = || {
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

        [block1, block2, block3]
    };

    let cfg = default_config();
    let window = Window {
        chr: "ref".to_string(),
        start: 1,
        end: 4,
    };

    let mut outputs: Vec<(Vec<String>, Vec<String>)> = Vec::new();
    for perm in perms {
        let mut all = build_blocks();
        let ordered: Vec<FasBlock> = perm
            .iter()
            .map(|&i| std::mem::replace(&mut all[i], make_block(vec![])))
            .collect();
        let blocks_per_input: Vec<Vec<FasBlock>> = ordered.into_iter().map(|b| vec![b]).collect();
        let merged = merge_window("ref", &window, &blocks_per_input, &cfg)
            .unwrap()
            .into_iter()
            .next()
            .expect("merge should succeed");
        let seqs: Vec<String> = merged
            .entries
            .iter()
            .map(|e| String::from_utf8(e.seq().to_vec()).unwrap())
            .collect();
        outputs.push((merged.names, seqs));
    }

    for (names, seqs) in &outputs[1..] {
        assert_eq!(names, &outputs[0].0);
        assert_eq!(seqs, &outputs[0].1);
    }
}

#[test]
fn merge_window_preserves_species_content() {
    // Two blocks with the same genome reference but different gap placements;
    // every species must keep all its ungapped bases after the merge.
    // MAF coordinate semantics: the reference range spans its non-gap bases
    // (4 bases here), so both blocks cover reference positions 1-4.
    let (ref_entry1, ref_name1, ref_header1) = make_entry("ref", 1, 4, "AC--GT");
    let (a_entry1, a_name1, a_header1) = make_entry("A", 1, 4, "AC--GT");
    let block1 = make_block(vec![
        (ref_entry1, ref_name1, ref_header1),
        (a_entry1, a_name1, a_header1),
    ]);

    let (ref_entry2, ref_name2, ref_header2) = make_entry("ref", 1, 4, "ACG-T-");
    let (b_entry2, b_name2, b_header2) = make_entry("B", 1, 6, "ACGTGT");
    let block2 = make_block(vec![
        (ref_entry2, ref_name2, ref_header2),
        (b_entry2, b_name2, b_header2),
    ]);

    let cfg = default_config();
    let window = Window {
        chr: "ref".to_string(),
        start: 1,
        end: 4,
    };

    let blocks_per_input = vec![vec![block1], vec![block2]];
    let merged = merge_window("ref", &window, &blocks_per_input, &cfg).unwrap();
    assert!(!merged.is_empty(), "merge should succeed");

    let total_ungapped = |name: &str| -> String {
        merged
            .iter()
            .flat_map(|block| {
                block
                    .entries
                    .iter()
                    .zip(block.names.iter())
                    .find(|(_, n)| n.as_str() == name)
                    .map(|(e, _)| {
                        e.seq()
                            .iter()
                            .filter(|&&b| b != b'-')
                            .copied()
                            .collect::<Vec<u8>>()
                    })
                    .unwrap_or_default()
            })
            .map(|b| b as char)
            .collect()
    };
    // Reference and A keep "ACGT" (4 bases); B keeps its 6 bases across the
    // merged overlap and the separately-emitted trailing insertion block.
    assert_eq!(total_ungapped("ref"), "ACGT");
    assert_eq!(total_ungapped("A"), "ACGT");
    assert_eq!(total_ungapped("B"), "ACGTGT");
}

#[test]
fn merge_window_union_of_multi_block_inputs() {
    // Real multi-pairwise style input: input 1 has two disjoint blocks
    // ([1-10], [21-30]) with species A; input 2 has one block [6-25] with
    // species B. The stream merge must keep the union of reference coverage
    // (1-30): front parts, merged overlaps and tails as separate blocks.
    let (ref1, n1, h1) = make_entry("ref", 1, 10, "ACGTACGTAC");
    let (a1, na1, ha1) = make_entry("A", 1, 10, "ACGTACGTAC");
    let block1a = make_block(vec![(ref1, n1, h1), (a1, na1, ha1)]);
    let (ref2, n2, h2) = make_entry("ref", 21, 30, "ACGTACGTAC");
    let (a2, na2, ha2) = make_entry("A", 21, 30, "ACGTACGTAC");
    let block1b = make_block(vec![(ref2, n2, h2), (a2, na2, ha2)]);
    // Same reference build: position p carries "ACGT"[(p-1) % 4], so block 2
    // (range 6-25) starts with CGTAC... and matches block 1 at positions 6-10.
    let (ref3, n3, h3) = make_entry("ref", 6, 25, "CGTACGTACGTACGTACGTA");
    let (b3, nb3, hb3) = make_entry("B", 6, 25, "CGTACGTACGTACGTACGTA");
    let block2 = make_block(vec![(ref3, n3, h3), (b3, nb3, hb3)]);

    let cfg = default_config();
    let window = Window {
        chr: "ref".to_string(),
        start: 1,
        end: 30,
    };
    let blocks_per_input = vec![vec![block1a, block1b], vec![block2]];
    let merged = merge_window("ref", &window, &blocks_per_input, &cfg).unwrap();

    // Union reference coverage 1-30 must be preserved across the output
    // blocks, and every species must keep its ungapped bases.
    let ref_cov: u64 = merged
        .iter()
        .map(|b| {
            let e = find_ref_entry(b, "ref").unwrap();
            (*e.range().end() as u64) - (*e.range().start() as u64) + 1
        })
        .sum();
    assert_eq!(ref_cov, 30, "union coverage 1-30, got {ref_cov}");

    let total_ungapped = |name: &str| -> usize {
        merged
            .iter()
            .map(|block| {
                block
                    .entries
                    .iter()
                    .zip(block.names.iter())
                    .find(|(_, n)| n.as_str() == name)
                    .map(|(e, _)| e.seq().iter().filter(|&&b| b != b'-').count())
                    .unwrap_or(0)
            })
            .sum()
    };
    assert_eq!(total_ungapped("A"), 20);
    assert_eq!(total_ungapped("B"), 20);
}

#[test]
fn merge_window_disjoint_blocks_both_kept() {
    // Two inputs whose blocks do not overlap at all: the stream merge emits
    // both as-is instead of dropping one.
    let (ref1, n1, h1) = make_entry("ref", 1, 10, "ACGTACGTAC");
    let (a1, na1, ha1) = make_entry("A", 1, 10, "ACGTACGTAC");
    let block1 = make_block(vec![(ref1, n1, h1), (a1, na1, ha1)]);
    let (ref2, n2, h2) = make_entry("ref", 20, 30, "ACGTACGTAC");
    let (b2, nb2, hb2) = make_entry("B", 20, 30, "ACGTACGTAC");
    let block2 = make_block(vec![(ref2, n2, h2), (b2, nb2, hb2)]);

    let cfg = default_config();
    let window = Window {
        chr: "ref".to_string(),
        start: 1,
        end: 30,
    };
    let blocks_per_input = vec![vec![block1], vec![block2]];
    let merged = merge_window("ref", &window, &blocks_per_input, &cfg).unwrap();
    assert_eq!(merged.len(), 2, "both disjoint blocks kept");
    let starts: Vec<i32> = merged
        .iter()
        .map(|b| *find_ref_entry(b, "ref").unwrap().range().start())
        .collect();
    assert_eq!(starts, vec![1, 20]);
}
