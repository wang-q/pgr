#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::io::Write;
use tempfile::Builder;

#[test]
fn command_nwk_cmp_single_file() {
    // Create a temporary Newick file with 2 different trees
    let mut file = Builder::new().suffix(".nwk").tempfile().unwrap();
    // Tree 1: ((A,B),(C,D)); -> Splits: {A,B} vs {C,D}
    // Tree 2: ((A,C),(B,D)); -> Splits: {A,C} vs {B,D}
    // RF distance should be 2.
    // Lengths are missing -> 0.0. WRF=0, KF=0.
    writeln!(file, "((A,B),(C,D));").unwrap();
    writeln!(file, "((A,C),(B,D));").unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "cmp", file.path().to_str().unwrap()])
        .run();

    // Expected output:
    // Tree1 Tree2 RF_Dist WRF_Dist KF_Dist
    // 1     1     0       0.000000 0.000000
    // 1     2     2       0.000000 0.000000
    // 2     1     2       0.000000 0.000000
    // 2     2     0       0.000000 0.000000

    assert!(stdout.contains("Tree1\tTree2\tRF_Dist\tWRF_Dist\tKF_Dist"));
    assert!(stdout.contains("1\t1\t0\t0\t0"));
    assert!(stdout.contains("1\t2\t2\t0\t0"));
    assert!(stdout.contains("2\t1\t2\t0\t0"));
    assert!(stdout.contains("2\t2\t0\t0\t0"));
}

#[test]
fn command_nwk_cmp_two_files() {
    let mut file1 = Builder::new().suffix(".nwk").tempfile().unwrap();
    writeln!(file1, "((A,B),(C,D));").unwrap(); // Tree 1

    let mut file2 = Builder::new().suffix(".nwk").tempfile().unwrap();
    writeln!(file2, "((A,B),(C,D));").unwrap(); // Tree 1 (Same)
    writeln!(file2, "((A,C),(B,D));").unwrap(); // Tree 2 (Diff, RF=2)

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "cmp",
            file1.path().to_str().unwrap(),
            file2.path().to_str().unwrap(),
        ])
        .run();

    // Expected:
    // T1(File1) vs T1(File2) -> 0
    // T1(File1) vs T2(File2) -> 2

    assert!(stdout.contains("1\t1\t0\t0\t0"));
    assert!(stdout.contains("1\t2\t2\t0\t0"));
}

#[test]
fn command_nwk_cmp_branch_lengths() {
    let mut file = Builder::new().suffix(".nwk").tempfile().unwrap();

    // T1: Same topology, lengths 0.2
    writeln!(file, "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);").unwrap();

    // T2: Same topology, one length 0.3 (Diff 0.1)
    writeln!(file, "((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.2);").unwrap();

    // T3: Diff topology, lengths 0.2
    writeln!(file, "((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);").unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "cmp", file.path().to_str().unwrap()])
        .run();

    // T1 vs T2: RF=0, WRF=0.1, KF=0.1
    // T1 vs T3: RF=2, WRF=0.8, KF=0.4

    // Check T1 vs T2
    // 1\t2\t0\t0.1\t0.1
    assert!(stdout.contains("1\t2\t0\t0.1\t0.1"));

    // Check T1 vs T3
    // 1\t3\t2\t0.8\t0.565685
    assert!(stdout.contains("1\t3\t2\t0.8\t0.565685"));
}
