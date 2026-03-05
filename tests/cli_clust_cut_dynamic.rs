use std::fs;
use tempfile::Builder;

#[path = "common/mod.rs"]
mod common;
use common::PgrCmd;

#[test]
fn test_dynamic_tree_cut_basic() {
    let temp = Builder::new().prefix("pgr_test_dynamic").tempdir().unwrap();
    let tree_file = temp.path().join("basic.nwk");
    let tree_content = "((A:0.1,B:0.1):0.8,(C:0.1,D:0.1):0.8);";
    fs::write(&tree_file, tree_content).unwrap();

    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "clust",
            "cut",
            tree_file.to_str().unwrap(),
            "--dynamic-tree",
            "2",
        ])
        .run();

    if !stderr.is_empty() {
        println!("STDERR: {}", stderr);
    }

    let lines: Vec<&str> = stdout.lines().collect();
    // Should have 2 clusters: {A,B} and {C,D}
    // Output format is tab separated members (with first being rep)

    let has_ab = lines.iter().any(|l| l.contains("A") && l.contains("B"));
    let has_cd = lines.iter().any(|l| l.contains("C") && l.contains("D"));

    assert!(has_ab, "Cluster {{A,B}} missing in output:\n{}", stdout);
    assert!(has_cd, "Cluster {{C,D}} missing in output:\n{}", stdout);
}

#[test]
fn test_dynamic_tree_cut_unassigned() {
    let temp = Builder::new()
        .prefix("pgr_test_dynamic_un")
        .tempdir()
        .unwrap();
    let tree_file = temp.path().join("unassigned.nwk");

    // Tree where min size is too large for leaves
    let tree_content = "((A:0.1,B:0.1):0.8,(C:0.1,D:0.1):0.8);";
    fs::write(&tree_file, tree_content).unwrap();

    // Min size 5. Total leaves 4.
    // Should result in unassigned (Cluster 0) or empty output if 0 is suppressed?
    // Our implementation currently outputs all clusters in the map.
    // Dynamic tree assigns 0 to unassigned nodes.
    // Partition.get_clusters() groups by value.
    // So we expect a cluster with ID 0 containing A,B,C,D (or however many are unassigned).
    // Or maybe multiple unassigned clusters? No, 0 is a single ID.

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "clust",
            "cut",
            tree_file.to_str().unwrap(),
            "--dynamic-tree",
            "5",
            "--format",
            "pair", // Easier to check IDs
        ])
        .run();

    // In pair format: Rep\tMember
    // If cluster ID is 0, it will be treated as a valid cluster by the writer code.
    // So we should see output.

    assert!(!stdout.is_empty());
    // Since all are unassigned (ID 0), they form one "cluster".
    // So we expect 4 lines, all belonging to the same representative.
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 4);
}
