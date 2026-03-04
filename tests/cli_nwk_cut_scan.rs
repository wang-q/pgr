use crate::common::*;
use std::fs;

mod common;

#[test]
fn test_scan_height() {
    // Tree: ((A:0.1,B:0.1):0.1,C:0.2);
    // Node heights (distance from leaves):
    // A, B, C: 0.0
    // (A,B): 0.1
    // Root ((A,B),C): 0.2
    let nwk = "((A:0.1,B:0.1):0.1,C:0.2);";
    let nwk_file = "tests/mat/scan_test.nwk";
    fs::write(nwk_file, nwk).expect("Failed to write nwk");

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "cut",
            nwk_file,
            "--height",
            "0",
            "--scan",
            "0,0.2,0.1",
        ])
        .run();

    let lines: Vec<&str> = stdout.lines().collect();
    // Header + 3 rows
    assert_eq!(lines.len(), 4, "Expected 4 lines output");
    assert_eq!(
        lines[0],
        "Threshold\tClusters\tSingletons\tNon-Singletons\tMaxSize"
    );

    // t=0
    let row0: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(row0[0], "0");
    assert_eq!(row0[1], "3"); // Clusters
    assert_eq!(row0[2], "3"); // Singletons
    assert_eq!(row0[3], "0"); // Non-Single
    assert_eq!(row0[4], "1"); // MaxSize

    // t=0.1
    let row1: Vec<&str> = lines[2].split('\t').collect();
    assert_eq!(row1[0], "0.1");
    assert_eq!(row1[1], "2");
    assert_eq!(row1[2], "1"); // {C}
    assert_eq!(row1[3], "1"); // {A,B}
    assert_eq!(row1[4], "2");

    // t=0.2
    let row2: Vec<&str> = lines[3].split('\t').collect();
    assert_eq!(row2[0], "0.2");
    assert_eq!(row2[1], "1");
    assert_eq!(row2[2], "0");
    assert_eq!(row2[3], "1"); // {A,B,C}
    assert_eq!(row2[4], "3");
}
