use crate::common::*;
use std::collections::HashSet;
use std::fs;

mod common;

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

fn parse_clusters(output: &str) -> Vec<HashSet<String>> {
    let mut clusters = Vec::new();
    for line in output.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let parts: HashSet<String> = line.split_whitespace().map(|s| s.to_string()).collect();
        clusters.push(parts);
    }
    clusters
}
