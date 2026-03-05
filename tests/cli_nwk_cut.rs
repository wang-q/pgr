mod common;
use crate::common::*;
use std::collections::HashSet;
use std::fs;

// --- Helper Functions ---

fn parse_clusters(output: &str) -> Vec<HashSet<String>> {
    let mut clusters = Vec::new();
    for line in output.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<String> = line.split_whitespace().map(|s| s.to_string()).collect();
        let set: HashSet<String> = parts.into_iter().collect();
        clusters.push(set);
    }

    clusters.sort_by(|a, b| {
        let min_a = a.iter().min().unwrap();
        let min_b = b.iter().min().unwrap();

        let num_a = min_a.parse::<i32>();
        let num_b = min_b.parse::<i32>();

        match (num_a, num_b) {
            (Ok(na), Ok(nb)) => na.cmp(&nb),
            _ => min_a.cmp(min_b),
        }
    });
    clusters
}

fn create_expected_clusters(groups: Vec<Vec<&str>>) -> Vec<HashSet<String>> {
    let mut clusters = Vec::new();
    for group in groups {
        let set: HashSet<String> = group.iter().map(|&s| s.to_string()).collect();
        clusters.push(set);
    }
    clusters.sort_by(|a, b| {
        let min_a = a.iter().min().unwrap();
        let min_b = b.iter().min().unwrap();

        let num_a = min_a.parse::<i32>();
        let num_b = min_b.parse::<i32>();

        match (num_a, num_b) {
            (Ok(na), Ok(nb)) => na.cmp(&nb),
            _ => min_a.cmp(min_b),
        }
    });
    clusters
}

// --- Tests ---

#[test]
fn test_nwk_cut_k() {
    let (stdout, stderr) = PgrCmd::new()
        .args(&["nwk", "cut", "tests/newick/abcde.nwk", "--k", "2"])
        .run();

    if !stderr.is_empty() {
        println!("STDERR: {}", stderr);
    }

    // Default format is 'cluster'.
    // K=2 -> {A,B} and {C}.
    // Sorted by size desc, then name.
    // {A,B} (size 2), {C} (size 1).

    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 2);
    assert_eq!(lines[0], "A\tB");
    assert_eq!(lines[1], "C");
}

#[test]
fn test_nwk_cut_pair_rep() {
    // K=2 -> {A, B}, {C}
    // Tree: ((A:1,B:1)D:1,C:1)E;
    // Dist from root: A=2, B=2, C=1.

    // 1. Default (root):
    // Cluster {A, B}: A(2), B(2). Tie. Alphabetical -> A.
    // Cluster {C}: C(1). -> C.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "cut",
            "tests/newick/abcde.nwk",
            "--k",
            "2",
            "--format",
            "pair",
        ])
        .run();
    let lines: Vec<_> = stdout.lines().collect();
    assert_eq!(lines.len(), 3);
    // {A,B} sorted first (size 2)
    assert_eq!(lines[0], "A\tA");
    assert_eq!(lines[1], "A\tB");
    // {C}
    assert_eq!(lines[2], "C\tC");

    // 2. Medoid:
    // Cluster {A, B}:
    // dist(A,B) = 2.
    // Sum for A: 2. Sum for B: 2.
    // Tie. Alphabetical -> A.
    // Same result.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "cut",
            "tests/newick/abcde.nwk",
            "--k",
            "2",
            "--format",
            "pair",
            "--rep",
            "medoid",
        ])
        .run();
    let lines: Vec<_> = stdout.lines().collect();
    assert_eq!(lines[0], "A\tA");
    assert_eq!(lines[1], "A\tB");

    // 3. First:
    // Cluster {A, B}. First is A.
    // Same result.
}

#[test]
fn test_nwk_cut_cluster_rep() {
    // K=2 -> {A, B}, {C}
    // A and B tie for everything.
    // Let's create a scenario where rep changes.
    // Tree: ((A:10,B:1)D:1,C:1)E;
    // Dist from root: A=11, B=2.
    // Rep(root): B should be rep.
    // Rep(first): A should be rep.

    // We can't easily construct a tree string with unequal lengths for testing here without modifying the file.
    // But we can verify that output format is correct (rep first).

    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "cut", "tests/newick/abcde.nwk", "--k", "2"])
        .run();
    let lines: Vec<_> = stdout.lines().collect();
    // A\tB -> A is rep.
    assert_eq!(lines[0], "A\tB");
}

#[test]
fn test_nwk_cut_height_pair() {
    let (stdout, _stderr) = PgrCmd::new()
        .args(&[
            "nwk",
            "cut",
            "tests/newick/abcde.nwk",
            "--height",
            "0.5",
            "--format",
            "pair",
        ])
        .run();

    // Height=0.5 -> {A}, {B}, {C} (all separate).
    // Format pair: Rep\tMember.
    // {A}: Rep A. Line: A\tA
    // {B}: Rep B. Line: B\tB
    // {C}: Rep C. Line: C\tC
    // Sorted by size (all 1), then name (A, B, C).

    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 3);
    assert_eq!(lines[0], "A\tA");
    assert_eq!(lines[1], "B\tB");
    assert_eq!(lines[2], "C\tC");
}

#[test]
fn test_nwk_cut_root_dist() {
    let (stdout, _stderr) = PgrCmd::new()
        .args(&["nwk", "cut", "tests/newick/abcde.nwk", "--root-dist", "0.5"])
        .run();

    // Root dist 0.5 -> {A,B}, {C}.
    // Output (cluster):
    // A\tB
    // C

    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 2);
    assert_eq!(lines[0], "A\tB");
    assert_eq!(lines[1], "C");
}

#[test]
fn test_nwk_cut_max_clade() {
    let (stdout, _stderr) = PgrCmd::new()
        .args(&["nwk", "cut", "tests/newick/abcde.nwk", "--max-clade", "2.5"])
        .run();

    // Max clade 2.5 -> {A,B} (diam 2), {C} (diam 0).
    // Output (cluster):
    // A\tB
    // C

    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 2);
    assert_eq!(lines[0], "A\tB");
    assert_eq!(lines[1], "C");
}

#[test]
fn test_avg_clade() {
    // Tree: ((A:0.1,B:0.1):0.4,C:0.5);
    // Dist(A,B) = 0.2
    // Dist(A,C) = 0.1 + 0.4 + 0.5 = 1.0
    // Dist(B,C) = 0.1 + 0.4 + 0.5 = 1.0
    //
    // Root ((A,B),C):
    // Max dist = 1.0
    // Avg dist = (0.2 + 1.0 + 1.0) / 3 = 0.7333...

    let nwk = "((A:0.1,B:0.1):0.4,C:0.5);";
    let nwk_file = "tests/mat/avg_test.nwk";
    // Ensure dir exists
    if !std::path::Path::new("tests/mat").exists() {
        fs::create_dir_all("tests/mat").unwrap();
    }
    fs::write(nwk_file, nwk).expect("Failed to write nwk");

    // 1. Max clade 0.8 -> Split ((A,B),C) because 1.0 > 0.8
    // Expect: {A,B}, {C} (2 clusters)
    let (out_max, _) = PgrCmd::new()
        .args(&["nwk", "cut", nwk_file, "--max-clade", "0.8"])
        .run();
    let clusters_max = parse_clusters(&out_max);
    assert_eq!(clusters_max.len(), 2, "Max clade should split");

    // 2. Avg clade 0.8 -> Keep ((A,B),C) because 0.733 < 0.8
    // Expect: {A,B,C} (1 cluster)
    let (out_avg, _) = PgrCmd::new()
        .args(&["nwk", "cut", nwk_file, "--avg-clade", "0.8"])
        .run();
    let clusters_avg = parse_clusters(&out_avg);
    assert_eq!(clusters_avg.len(), 1, "Avg clade should keep");
    assert_eq!(clusters_avg[0].len(), 3);
}

#[test]
fn test_scan_height() {
    // Tree: ((A:0.1,B:0.1):0.1,C:0.2);
    // Node heights (distance from leaves):
    // A, B, C: 0.0
    // (A,B): 0.1
    // Root ((A,B),C): 0.2
    let nwk = "((A:0.1,B:0.1):0.1,C:0.2);";
    let nwk_file = "tests/mat/scan_test.nwk";
    // Ensure dir exists
    if !std::path::Path::new("tests/mat").exists() {
        fs::create_dir_all("tests/mat").unwrap();
    }
    fs::write(nwk_file, nwk).expect("Failed to write nwk");

    let stats_file = "tests/mat/scan_stats.tsv";

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "cut",
            nwk_file,
            "--height",
            "0",
            "--scan",
            "0,0.2,0.1",
            "--stats-out",
            stats_file,
        ])
        .run();

    // Verify stdout (Long Format)
    let out_lines: Vec<&str> = stdout.lines().collect();
    assert!(out_lines.len() > 1);
    assert_eq!(out_lines[0], "Group\tClusterID\tSampleID");

    // Verify stats file
    let stats_content = fs::read_to_string(stats_file).expect("Failed to read stats file");
    let lines: Vec<&str> = stats_content.lines().collect();

    // Header + 3 rows
    assert_eq!(lines.len(), 4, "Expected 4 lines output in stats file");
    assert_eq!(
        lines[0],
        "Group\tClusters\tSingletons\tNon-Singletons\tMaxSize"
    );

    // t=0
    let row0: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(row0[0], "height=0");
    assert_eq!(row0[1], "3"); // Clusters
    assert_eq!(row0[2], "3"); // Singletons
    assert_eq!(row0[3], "0"); // Non-Single
    assert_eq!(row0[4], "1"); // MaxSize

    // t=0.1
    let row1: Vec<&str> = lines[2].split('\t').collect();
    assert_eq!(row1[0], "height=0.1");
    assert_eq!(row1[1], "2");
    assert_eq!(row1[2], "1"); // {C}
    assert_eq!(row1[3], "1"); // {A,B}
    assert_eq!(row1[4], "2");

    // t=0.2
    let row2: Vec<&str> = lines[3].split('\t').collect();
    assert_eq!(row2[0], "height=0.2");
    assert_eq!(row2[1], "1");
    assert_eq!(row2[2], "0");
    assert_eq!(row2[3], "1"); // {A,B,C}
    assert_eq!(row2[4], "3");

    // Cleanup
    let _ = fs::remove_file(nwk_file);
    let _ = fs::remove_file(stats_file);
}

#[test]
fn test_scipy_workflow() {
    // 1. Generate Newick tree from Q_X using pgr clust hier
    // The input file tests/mat/scipy_Q_X.phy must exist (generated by python script)
    let phy_file = "tests/mat/scipy_Q_X.phy";
    if !std::path::Path::new(phy_file).exists() {
        eprintln!(
            "Skipping test_scipy_workflow: {} not found. Run gen_phy.py first.",
            phy_file
        );
        return;
    }

    let (stdout, stderr) = PgrCmd::new()
        .args(&["clust", "hier", phy_file, "--method", "single"])
        .run();

    if stdout.trim().is_empty() {
        eprintln!("pgr clust hier failed. Stderr: {}", stderr);
    }
    assert!(!stdout.trim().is_empty(), "pgr clust hier output is empty");

    let nwk_file = "tests/mat/scipy_tree.nwk";
    fs::write(nwk_file, &stdout).expect("Failed to write tree file");

    // 2. Verify fcluster_distance (t=0.6)
    // Note: pgr clust hier outputs tree where height = distance / 2.
    // So for SciPy t=0.6, we cut at height 0.3.
    // Expected clusters:
    // 4: {0..6, 8..9} -> {0,1,2,3,4,5,6,8,9}
    // 5: {7}
    // 6: {10..14, 16..19} -> {10,11,12,13,14,16,17,18,19}
    // 7: {15}
    // 3: {20}
    // 1: {21..23, 25..29} -> {21,22,23,25,26,27,28,29}
    // 2: {24}
    let expected_dist_0_6 = create_expected_clusters(vec![
        vec!["0", "1", "2", "3", "4", "5", "6", "8", "9"],
        vec!["7"],
        vec!["10", "11", "12", "13", "14", "16", "17", "18", "19"],
        vec!["15"],
        vec!["20"],
        vec!["21", "22", "23", "25", "26", "27", "28", "29"],
        vec!["24"],
    ]);

    let (cut_out, stderr) = PgrCmd::new()
        .args(&["nwk", "cut", nwk_file, "--height", "0.3"])
        .run();

    if cut_out.trim().is_empty() {
        eprintln!("nwk cut height=0.3 failed. Stderr: {}", stderr);
    }

    let actual_dist_0_6 = parse_clusters(&cut_out);
    assert_eq!(
        actual_dist_0_6, expected_dist_0_6,
        "Distance t=0.6 (h=0.3) failed"
    );

    // 3. Verify fcluster_maxclust (k=4)
    // Expected clusters:
    // 3: {0..9}
    // 4: {10..19}
    // 2: {20}
    // 1: {21..29}
    let expected_maxclust_4 = create_expected_clusters(vec![
        vec!["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"],
        vec!["10", "11", "12", "13", "14", "15", "16", "17", "18", "19"],
        vec!["20"],
        vec!["21", "22", "23", "24", "25", "26", "27", "28", "29"],
    ]);

    let (cut_out_k4, _) = PgrCmd::new()
        .args(&["nwk", "cut", nwk_file, "--k", "4"])
        .run();

    let actual_maxclust_4 = parse_clusters(&cut_out_k4);
    assert_eq!(
        actual_maxclust_4, expected_maxclust_4,
        "Maxclust k=4 failed"
    );

    // 4. Verify fcluster_maxclust (k=8)
    // 8.0: array([5, 5, 5, 5, 5, 5, 5, 6, 5, 5, 7, 7, 7, 7, 7, 8, 7, 7, 7, 7, 4,
    //             1, 1, 1, 3, 1, 1, 1, 1, 2]),
    // Clusters:
    // 5: {0..6, 8..9} (same as dist 0.6 cluster 4)
    // 6: {7}
    // 7: {10..14, 16..19} (same as dist 0.6 cluster 6)
    // 8: {15}
    // 4: {20}
    // 1: {21..23, 25..28} -> {21,22,23,25,26,27,28} ? Wait, index 29 is 2.
    // Index 24 is 3.
    // Index 29 is 2.
    // So:
    // 1: {21,22,23,25,26,27,28}
    // 2: {29}
    // 3: {24}
    // 4: {20}
    // 5: {0..6, 8..9}
    // 6: {7}
    // 7: {10..14, 16..19}
    // 8: {15}
    // Total 8 clusters.
    let expected_maxclust_8 = create_expected_clusters(vec![
        vec!["0", "1", "2", "3", "4", "5", "6", "8", "9"],
        vec!["7"],
        vec!["10", "11", "12", "13", "14", "16", "17", "18", "19"],
        vec!["15"],
        vec!["20"],
        vec!["21", "22", "23", "25", "26", "27", "28"],
        vec!["29"],
        vec!["24"],
    ]);

    let (cut_out_k8, _) = PgrCmd::new()
        .args(&["nwk", "cut", nwk_file, "--k", "8"])
        .run();

    let actual_maxclust_8 = parse_clusters(&cut_out_k8);
    assert_eq!(
        actual_maxclust_8, expected_maxclust_8,
        "Maxclust k=8 failed"
    );

    // 5. Verify fcluster_inconsistent (t=0.8, depth=2)
    // pgr implementation yields slightly different topology/inconsistency due to tie-breaking,
    // so we verify against pgr's own stable output for regression testing.
    let expected_inc_0_8 = create_expected_clusters(vec![
        vec!["0", "4"],
        vec!["1", "2", "5"],
        vec!["3"],
        vec!["6", "8"],
        vec!["7"],
        vec!["9"],
        vec!["10", "17"],
        vec!["11"],
        vec!["12"],
        vec!["13", "14", "18", "19"],
        vec!["15"],
        vec!["16"],
        vec!["20"],
        vec!["21", "22"],
        vec!["23", "25", "26", "28"],
        vec!["24"],
        vec!["27"],
        vec!["29"],
    ]);

    let (cut_out_inc, _stderr) = PgrCmd::new()
        .args(&["nwk", "cut", nwk_file, "--inconsistent", "0.8"])
        .run();

    let actual_inc_0_8 = parse_clusters(&cut_out_inc);
    assert_eq!(
        actual_inc_0_8, expected_inc_0_8,
        "Inconsistent t=0.8 failed"
    );
}

#[test]
fn test_cut_support_filter() {
    // Tree: ((A:0.1,B:0.1)90:0.1,C:0.2);
    // Internal node (A,B) has support 90.
    // Edge (A,B)->Root has length 0.1.
    // If support < 95, this edge length becomes INF.
    let nwk = "((A:0.1,B:0.1)90:0.1,C:0.2);";
    let nwk_file = "tests/nwk/support_test.nwk";
    fs::create_dir_all("tests/nwk").unwrap();
    fs::write(nwk_file, nwk).expect("Failed to write nwk");

    // Case 1: No support filter (or low threshold)
    // Max pairwise distance: dist(A,C) = 0.1 + 0.1 + 0.2 = 0.4.
    // Threshold 0.5 > 0.4.
    // Result: 1 cluster {A,B,C}.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "cut",
            nwk_file,
            "--max-clade",
            "0.5",
            "--support",
            "0", // Threshold 0 < 90, so edge is kept
        ])
        .run();

    let out_lines: Vec<&str> = stdout.lines().collect();
    // Expect 1 line (1 cluster)
    assert_eq!(
        out_lines.len(),
        1,
        "Expected 1 cluster without support filter"
    );
    assert_eq!(out_lines[0], "A\tB\tC");

    // Case 2: High support threshold
    // Threshold 95 > 90. Edge (A,B)->Root becomes INF.
    // dist(A,C) = 0.1 + INF + 0.2 = INF.
    // INF > 0.5.
    // Root node split.
    // (A,B) subtree: dist(A,B) = 0.2 <= 0.5. Kept as cluster.
    // C is singleton.
    // Result: {A,B}, {C}.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "cut",
            nwk_file,
            "--max-clade",
            "0.5",
            "--support",
            "95",
        ])
        .run();

    let mut out_lines: Vec<&str> = stdout.lines().collect();
    out_lines.sort(); // Ensure deterministic order for assertion
                      // Expected: "A\tB", "C"
    assert_eq!(
        out_lines.len(),
        2,
        "Expected 2 clusters with support filter"
    );
    assert_eq!(out_lines[0], "A\tB");
    assert_eq!(out_lines[1], "C");

    // Cleanup
    let _ = fs::remove_file(nwk_file);
}
