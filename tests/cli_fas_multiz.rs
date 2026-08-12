#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;
use tempfile::TempDir;

#[test]
fn command_fas_multiz_merges_inputs() {
    let tempdir = TempDir::new().unwrap();
    let out_path = tempdir.path().join("merged.fas");
    let out_str = out_path.to_str().unwrap();

    PgrCmd::new()
        .args(&[
            "fas",
            "multiz",
            "-r",
            "S288c",
            "tests/fas/S288cvsRM11_1a.slice.fas",
            "tests/fas/S288cvsYJM789.slice.fas",
            "tests/fas/S288cvsSpar.slice.fas",
            "-o",
            out_str,
        ])
        .assert()
        .success();

    assert!(out_path.is_file());
    let content = fs::read_to_string(out_path).unwrap();
    assert!(content.lines().count() > 0);

    tempdir.close().unwrap();
}

#[test]
fn command_fas_multiz_merges_conflicting_refs() {
    let tempdir = TempDir::new().unwrap();

    // Two pairwise blocks sharing the S288c reference, but the reference
    // sequence differs at one SNP position in the second input. Without the
    // crossover merge this window would be dropped; with it, the merge
    // succeeds by splicing the better-matching side.
    let seq_a = "ACGT".repeat(20);
    let seq_b = format!("{}{}{}", &seq_a[..16], "T", &seq_a[17..]);
    let a_path = tempdir.path().join("a.fas");
    let b_path = tempdir.path().join("b.fas");
    fs::write(
        &a_path,
        format!(
            ">S288c.I(+):100-179\n{}\n>RM11_1a.I(+):500-579\n{}\n",
            seq_a, seq_a
        ),
    )
    .unwrap();
    fs::write(
        &b_path,
        format!(
            ">S288c.I(+):100-179\n{}\n>RM11_1a.I(+):500-579\n{}\n",
            seq_b, seq_a
        ),
    )
    .unwrap();

    let out_path = tempdir.path().join("merged.fas");
    let out_str = out_path.to_str().unwrap();
    let a_str = a_path.to_str().unwrap();
    let b_str = b_path.to_str().unwrap();

    PgrCmd::new()
        .args(&["fas", "multiz", "-r", "S288c", a_str, b_str, "-o", out_str])
        .assert()
        .success();

    let content = fs::read_to_string(&out_path).unwrap();
    // Merged block must contain both species.
    assert!(content.contains(">S288c.I(+):100-179"));
    assert!(content.contains(">RM11_1a.I(+):500-579"));

    tempdir.close().unwrap();
}
