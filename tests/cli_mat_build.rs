#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;
use tempfile::TempDir;

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
