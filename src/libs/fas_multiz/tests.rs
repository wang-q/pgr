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
    use intspan::Range;
    use std::fs::File;
    use std::io::Write;

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
        assert_eq!(
            names,
            vec!["A".to_string(), "B".to_string(), "ref".to_string()]
        );
    }
}

#[test]
fn merge_fas_files_auto_windows_matches_explicit() {
    use intspan::Range;
    use std::fs::File;
    use std::io::Write;

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
        vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "ref".to_string()
        ]
    );
}
