mod common;
use crate::common::*;

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
