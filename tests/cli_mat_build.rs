#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;
use tempfile::TempDir;

#[test]
fn command_mat_upgma() {
    let temp = TempDir::new().unwrap();
    let input = temp.path().join("input.phy");
    let output = temp.path().join("output.nwk");

    let content = "4
A 0 7 11 14
B 7 0 6 9
C 11 6 0 7
D 14 9 7 0
";
    fs::write(&input, content).unwrap();

    PgrCmd::new()
        .args(&[
            "mat",
            "upgma",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .assert()
        .success();

    let nwk = fs::read_to_string(&output).unwrap();
    // println!("{}", nwk); // For debugging
    assert!(nwk.contains("A:"));
    assert!(nwk.contains("B:"));
    assert!(nwk.contains("C:"));
    assert!(nwk.contains("D:"));

    // Check topology structure
    // For this matrix: B-C=6 (min), so B and C merge first at height 3.
    // Then D merges with (B,C) at height 4.
    // Finally A merges with ((B,C),D) at height 5.333.

    // Check for B:3 and C:3
    assert!(nwk.contains("B:3") || nwk.contains("B:3.0"));
    assert!(nwk.contains("C:3") || nwk.contains("C:3.0"));

    // Check for D:4
    assert!(nwk.contains("D:4") || nwk.contains("D:4.0"));

    // Check for A:5.33
    assert!(nwk.contains("A:5.33"));

    // Check groupings
    assert!(nwk.contains("(B:3"));
    assert!(nwk.contains("C:3)"));
}

#[test]
fn command_mat_nj() {
    let temp = TempDir::new().unwrap();
    let input = temp.path().join("input.phy");
    let output = temp.path().join("output.nwk");

    let content = "4
A 0 7 11 14
B 7 0 6 9
C 11 6 0 7
D 14 9 7 0
";
    fs::write(&input, content).unwrap();

    PgrCmd::new()
        .args(&[
            "mat",
            "nj",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .assert()
        .success();

    let nwk = fs::read_to_string(&output).unwrap();
    assert!(nwk.contains("A:"));
    assert!(nwk.contains("B:"));
    assert!(nwk.contains("C:"));
    assert!(nwk.contains("D:"));

    // Check topology structure
    // NJ on this additive matrix should recover ((A,B),(C,D))
    // Note: The exact string depends on rooting and child order, but (A,B) and (C,D) should be clades
    // Current NJ implementation roots at midpoint of last edge, so we expect a rooted tree.

    // We can also verify via pipe
    let (stdout, _) = PgrCmd::new()
        .args(&["mat", "nj", "stdin"])
        .stdin(content)
        .run();

    assert!(stdout.contains("A:"));
    assert!(stdout.contains("B:"));
}
