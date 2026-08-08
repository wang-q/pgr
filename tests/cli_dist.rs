#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::path::PathBuf;

/// Return the absolute path to a fixture in `tests/dist/input`.
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/dist/input")
        .join(name)
}

#[test]
fn command_dist_hv() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "hv",
            fixture("seq.fa").to_str().unwrap(),
            "-k",
            "7",
            "-w",
            "1",
        ])
        .run();

    assert!(stdout.lines().count() >= 1);
    assert!(stdout.contains(fixture("seq.fa").to_str().unwrap()));
}

#[test]
fn command_dist_hv_pair() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "hv",
            fixture("seq.fa").to_str().unwrap(),
            fixture("seq.fa").to_str().unwrap(), // Compare file against itself
        ])
        .run();

    assert!(stdout.contains(fixture("seq.fa").to_str().unwrap()));
    // Similarity should be 1.0 / Distance 0.0
    // The output format: <file1> <file2> ... <mash_dist> ...
}

#[test]
fn command_dist_hv_syncmer() {
    // Closed syncmers projected onto hypervectors; syng DNA defaults (smer=8,
    // window=55) applied automatically.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "hv",
            fixture("genome1.fa").to_str().unwrap(),
            fixture("genome1.fa").to_str().unwrap(),
            "--sampler",
            "syncmer",
        ])
        .run();

    // Single pair of files -> one line.
    assert_eq!(stdout.lines().count(), 1);
    // Self-comparison: identical merged syncmer HV -> distance 0, jaccard 1.
    assert!(stdout.contains("0.0000"));
    assert!(stdout.contains("1.0000"));
}

#[test]
fn command_dist_hv_syncmer_defaults() {
    // seq.fa (4 short periodic sequences) merged to one HV; DNA syncmer
    // defaults smer=8/window=55 produce 2 unique syncmers -> card=2.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "hv",
            fixture("seq.fa").to_str().unwrap(),
            "--sampler",
            "syncmer",
        ])
        .run();

    assert!(stdout.contains("\t2\t2\t2\t2\t0.0000\t1.0000\t1.0000"));
}

#[test]
fn command_dist_hv_syncmer_explicit_params() {
    // Explicit -k/-w change the sampled syncmer set (3 vs 2 unique).
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "hv",
            fixture("seq.fa").to_str().unwrap(),
            "--sampler",
            "syncmer",
            "-k",
            "5",
            "-w",
            "3",
        ])
        .run();

    assert!(stdout.contains("\t3\t3\t3\t3\t0.0000\t1.0000\t1.0000"));
}

#[test]
fn command_dist_hv_syncmer_protein_defaults() {
    // --protein applies protein syncmer defaults (smer=7, window=5), a
    // different sketch than the DNA defaults (1 vs 2 unique syncmers).
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "hv",
            fixture("seq.fa").to_str().unwrap(),
            "--sampler",
            "syncmer",
            "--protein",
        ])
        .run();

    assert!(stdout.contains("\t1\t1\t1\t1\t0.0000\t1.0000\t1.0000"));
}

#[test]
fn command_dist_hv_protein_rejects_mod_hasher() {
    // --hasher mod is DNA-only; with --protein it must error.
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "dist",
            "hv",
            fixture("seq.fa").to_str().unwrap(),
            "--hasher",
            "mod",
            "--protein",
        ])
        .run_fail();

    assert!(stderr.contains("--hasher mod is DNA-only"));
}

#[test]
fn command_dist_hv_files_self_and_dim_mismatch() {
    let temp = tempfile::TempDir::new().unwrap();
    let fa = temp.path().join("g.fa");
    std::fs::write(&fa, format!(">g\n{}\n", "ACGT".repeat(50))).unwrap();
    let idx = temp.path().join("g.pgi");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "pgi",
            "build",
            fa.to_str().unwrap(),
            "-o",
            idx.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("wrote"));

    let hv1 = temp.path().join("g.hv");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "pgi",
            "to-hv",
            idx.to_str().unwrap(),
            "-o",
            hv1.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("hypervector"));

    // Self-comparison of a .hv file: mash 0, jaccard 1.
    let (stdout, _) = PgrCmd::new()
        .args(&["dist", "hv", hv1.to_str().unwrap(), hv1.to_str().unwrap()])
        .run();
    let fields: Vec<&str> = stdout.split_whitespace().collect();
    assert_eq!(fields.len(), 9, "bad hv output: {stdout}");
    assert_eq!(fields[0], "g");
    assert_eq!(fields[6], "0.0000", "self mash must be 0: {stdout}");
    assert_eq!(fields[7], "1.0000", "self jaccard must be 1: {stdout}");

    // Different dimensions must fail loudly.
    let hv2 = temp.path().join("g512.hv");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "pgi",
            "to-hv",
            idx.to_str().unwrap(),
            "-o",
            hv2.to_str().unwrap(),
            "--dim",
            "512",
        ])
        .run();
    assert!(stderr.contains("hypervector"));
    let (_, stderr) = PgrCmd::new()
        .args(&["dist", "hv", hv1.to_str().unwrap(), hv2.to_str().unwrap()])
        .run_fail();
    assert!(
        stderr.contains("dimension mismatch"),
        "expected dim mismatch: {stderr}"
    );
}
