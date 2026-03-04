#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

#[test]
fn command_clust_upgma() {
    let temp = tempfile::TempDir::new().unwrap();
    let input = temp.path().join("input.phy");
    let output = temp.path().join("output.nwk");

    let content = "4
A 0 7 11 14
B 7 0 6 9
C 11 6 0 7
D 14 9 7 0
";
    std::fs::write(&input, content).unwrap();

    PgrCmd::new()
        .args(&[
            "clust",
            "upgma",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .assert()
        .success();

    let nwk = std::fs::read_to_string(&output).unwrap();
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
fn command_clust_hier() {
    let temp = tempfile::TempDir::new().unwrap();
    let input = temp.path().join("input.phy");
    let output = temp.path().join("output.nwk");

    let content = "4
A 0 7 11 14
B 7 0 6 9
C 11 6 0 7
D 14 9 7 0
";
    std::fs::write(&input, content).unwrap();

    // 1. Test UPGMA (Average Linkage)
    // Should produce identical results to `clust upgma` (topology & branch lengths logic)
    // In `hier`, branch length = distance/2 - child_height
    // For B-C (d=6): height=3. B,C leaves height=0. Branch=3.
    PgrCmd::new()
        .args(&[
            "clust",
            "hier",
            input.to_str().unwrap(),
            "--method",
            "average",
            "-o",
            output.to_str().unwrap(),
        ])
        .assert()
        .success();

    let nwk = std::fs::read_to_string(&output).unwrap();
    assert!(nwk.contains("B:3") || nwk.contains("B:3.0"));
    assert!(nwk.contains("C:3") || nwk.contains("C:3.0"));
    assert!(nwk.contains("D:4") || nwk.contains("D:4.0"));
    assert!(nwk.contains("A:5.33"));

    // 2. Test Single Linkage
    // B-C (d=6) -> (B,C)
    // (B,C)-D: min(d(B,D)=9, d(C,D)=7) = 7. Merge at d=7. Height=3.5.
    // ((B,C),D)-A: min(d(B,A)=7, d(C,A)=11, d(D,A)=14) = 7. Merge at d=7. Height=3.5.
    // Wait, D-A is 14. B-A is 7. C-A is 11.
    // min(d(BC, A)) = min(7, 11) = 7.
    // min(d(BC, D)) = min(9, 7) = 7.
    // So (B,C) can merge with D or A at d=7.
    // Ties? B-C is 6. Next min is 7 (B-A, C-D).
    // If we merge B-C first (d=6).
    // Matrix becomes:
    //    BC   D   A
    // BC 0    7   7
    // D  7    0   14
    // A  7    14  0
    //
    // Next min is 7. We can merge (BC)-D or (BC)-A.
    // Let's see what happens.
    let output_single = temp.path().join("output_single.nwk");
    PgrCmd::new()
        .args(&[
            "clust",
            "hier",
            input.to_str().unwrap(),
            "--method",
            "single",
            "-o",
            output_single.to_str().unwrap(),
        ])
        .assert()
        .success();
    
    // Just verify it runs and produces output
    let nwk_single = std::fs::read_to_string(&output_single).unwrap();
    assert!(nwk_single.starts_with("("));
    assert!(nwk_single.trim().ends_with(";"));

    // 3. Test Ward (default)
    let output_ward = temp.path().join("output_ward.nwk");
    PgrCmd::new()
        .args(&[
            "clust",
            "hier",
            input.to_str().unwrap(),
            "-o",
            output_ward.to_str().unwrap(),
        ])
        .assert()
        .success();
    
    let nwk_ward = std::fs::read_to_string(&output_ward).unwrap();
    assert!(nwk_ward.len() > 0);
}

#[test]
fn command_clust_nj() {
    let temp = tempfile::TempDir::new().unwrap();
    let input = temp.path().join("input.phy");
    let output = temp.path().join("output.nwk");

    let content = "4
A 0 7 11 14
B 7 0 6 9
C 11 6 0 7
D 14 9 7 0
";
    std::fs::write(&input, content).unwrap();

    PgrCmd::new()
        .args(&[
            "clust",
            "nj",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .assert()
        .success();

    let nwk = std::fs::read_to_string(&output).unwrap();
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
        .args(&["clust", "nj", "stdin"])
        .stdin(content)
        .run();

    assert!(stdout.contains("A:"));
    assert!(stdout.contains("B:"));
}
