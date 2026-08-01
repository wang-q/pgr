#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;

/// Write a FASTA with one sequence `seq` named `name` into `dir/name.fa`.
fn write_fa(dir: &std::path::Path, name: &str, seq: &str) -> String {
    let path = dir.join(format!("{name}.fa"));
    fs::write(&path, format!(">{name}\n{seq}\n")).unwrap();
    path.to_string_lossy().to_string()
}

#[test]
fn command_pgi_build_stat() {
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_fa(temp.path(), "g", &"ACGT".repeat(100));
    let out = temp.path().join("g.pgi");
    let (_, stderr) = PgrCmd::new()
        .args(&["pgi", "build", &fa, "-o", out.to_str().unwrap()])
        .run();
    assert!(stderr.contains("wrote"), "missing build log: {stderr}");
    assert!(out.exists());

    let (stdout, _) = PgrCmd::new()
        .args(&["pgi", "stat", out.to_str().unwrap()])
        .run();
    assert!(stdout.contains("K-mer size: 40"), "got {stdout}");
    assert!(stdout.contains("Syncmer: 8/5"), "got {stdout}");
    assert!(stdout.contains("Contigs: 1"), "got {stdout}");
    assert!(stdout.contains("Unique k-mers:"), "got {stdout}");
}

#[test]
fn command_pgi_build_from_2bit() {
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_fa(temp.path(), "g", &"ACGT".repeat(100));
    let tb = temp.path().join("g.2bit");
    let (_, stderr) = PgrCmd::new()
        .args(&["fa", "to-2bit", &fa, "-o", tb.to_str().unwrap()])
        .run();
    assert!(!stderr.contains("error:"), "2bit build failed: {stderr}");

    let out = temp.path().join("g.pgi");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "pgi",
            "build",
            tb.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(
        stderr.contains("wrote"),
        "pgi build from 2bit failed: {stderr}"
    );
    assert!(out.exists());
}

#[test]
fn command_pgi_dist_identical_and_disjoint() {
    let temp = tempfile::TempDir::new().unwrap();
    let fa1 = write_fa(temp.path(), "a", &"ACGT".repeat(100));
    let fa2 = write_fa(temp.path(), "a2", &"ACGT".repeat(100));
    let fa3 = write_fa(temp.path(), "b", &"TTTT".repeat(100));
    let idx1 = temp.path().join("a.pgi");
    let idx2 = temp.path().join("a2.pgi");
    let idx3 = temp.path().join("b.pgi");
    for (fa, out) in [(&fa1, &idx1), (&fa2, &idx2), (&fa3, &idx3)] {
        let (_, stderr) = PgrCmd::new()
            .args(&["pgi", "build", fa, "-o", out.to_str().unwrap()])
            .run();
        assert!(stderr.contains("wrote"), "build failed: {stderr}");
    }

    // Identical sequences -> Jaccard 1, Mash 0.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "pgi",
            idx1.to_str().unwrap(),
            idx2.to_str().unwrap(),
        ])
        .run();
    let fields: Vec<&str> = stdout.trim().split('\t').collect();
    assert_eq!(fields.len(), 9, "unexpected output: {stdout}");
    assert_eq!(
        fields[4], fields[5],
        "identical indexes should have inter == union"
    );
    assert_eq!(fields[7], "1.0000", "jaccard should be 1: {stdout}");

    // Disjoint sequences -> Jaccard 0.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "pgi",
            idx1.to_str().unwrap(),
            idx3.to_str().unwrap(),
        ])
        .run();
    let fields: Vec<&str> = stdout.trim().split('\t').collect();
    assert_eq!(fields[4], "0", "inter should be 0: {stdout}");
    assert_eq!(fields[7], "0.0000", "jaccard should be 0: {stdout}");
}

#[test]
fn command_pgi_dist_param_mismatch_fails() {
    let temp = tempfile::TempDir::new().unwrap();
    let fa1 = write_fa(temp.path(), "a", &"ACGT".repeat(100));
    let fa2 = write_fa(temp.path(), "b", &"ACGT".repeat(100));
    let idx1 = temp.path().join("a.pgi");
    let idx2 = temp.path().join("b.pgi");
    let (_, stderr) = PgrCmd::new()
        .args(&["pgi", "build", &fa1, "-o", idx1.to_str().unwrap()])
        .run();
    assert!(stderr.contains("wrote"));
    // Different k-mer size.
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "pgi",
            "build",
            &fa2,
            "-o",
            idx2.to_str().unwrap(),
            "--kmer",
            "20",
        ])
        .run();
    assert!(stderr.contains("wrote"));

    let (_, stderr) = PgrCmd::new()
        .args(&[
            "dist",
            "pgi",
            idx1.to_str().unwrap(),
            idx2.to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("k-mer size mismatch"),
        "expected mismatch error: {stderr}"
    );
}

#[test]
fn command_pgi_to_hv() {
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_fa(temp.path(), "g", &"ACGT".repeat(100));
    let idx = temp.path().join("g.pgi");
    let (_, stderr) = PgrCmd::new()
        .args(&["pgi", "build", &fa, "-o", idx.to_str().unwrap()])
        .run();
    assert!(stderr.contains("wrote"));

    let hv = temp.path().join("g.hv");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "pgi",
            "to-hv",
            idx.to_str().unwrap(),
            "-o",
            hv.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("hypervector"), "to-hv failed: {stderr}");
    let bytes = fs::read(&hv).unwrap();
    assert_eq!(&bytes[0..4], b"PGV1", "bad hv magic");
    // 4 magic + 4 ver + 4 k + 4 dim + 4 name_len + name(1) + 1024*4 hv bytes
    assert_eq!(bytes.len(), 4 + 4 + 4 + 4 + 4 + 1 + 1024 * 4);
}

#[test]
fn command_pgi_build_invalid_k_fails() {
    let temp = tempfile::TempDir::new().unwrap();
    let fa = write_fa(temp.path(), "g", &"ACGT".repeat(100));
    let out = temp.path().join("g.pgi");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "pgi",
            "build",
            &fa,
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "65",
        ])
        .run_fail();
    assert!(stderr.contains("k must be in 1..=64"), "got {stderr}");
}
