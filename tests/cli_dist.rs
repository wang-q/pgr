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
            fixture("genome2.fa").to_str().unwrap(),
            "--sampler",
            "syncmer",
        ])
        .run();

    // Single pair of files -> one line.
    assert_eq!(stdout.lines().count(), 1);
    assert!(stdout.contains("0.0000"));
    assert!(stdout.contains("1.0000"));
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

#[test]
fn command_dist_seq() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "seq",
            fixture("seq.fa").to_str().unwrap(),
            "-k",
            "7",
            "-w",
            "1",
            "--zero",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 16);
    assert!(stdout.contains("seqA\tseqB\t0.0168\t0.8000\t1.0000"));
}

#[test]
fn command_dist_seq_sim() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "seq",
            fixture("seq.fa").to_str().unwrap(),
            "-k",
            "7",
            "-w",
            "1",
            "--zero",
            "--sim",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 16);
    // Mash dist 0.0168 -> Sim 1 - 0.0168 = 0.9832
    assert!(stdout.contains("seqA\tseqB\t0.9832\t0.8000\t1.0000"));
}

#[test]
fn command_dist_seq_genome() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "seq",
            fixture("genome1.fa").to_str().unwrap(),
            fixture("genome2.fa").to_str().unwrap(),
            "-k",
            "21",
            "-w",
            "5",
            "--hasher",
            "mod",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 2);
    assert!(stdout.contains("chrA\tchrA\t0.0000\t1.0000\t1.0000"));
    assert!(stdout.contains("chrB\tchrA\t0.0597\t0.1667\t0.1667"));
}

#[test]
fn command_dist_seq_merge() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "seq",
            fixture("seq.fa").to_str().unwrap(),
            "-k",
            "7",
            "-w",
            "1",
            "--merge",
            "--hasher",
            "murmur",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 1);
    assert!(stdout.contains(&format!(
        "{}\t{}\t9\t9\t9\t9\t0.0000\t1.0000\t1.0000",
        fixture("seq.fa").to_str().unwrap(),
        fixture("seq.fa").to_str().unwrap()
    )));
}

#[test]
fn command_dist_seq_syncmer() {
    // Closed syncmers with syng DNA defaults (smer=8, window=55) applied
    // automatically when --sampler syncmer is used without explicit -k/-w.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "seq",
            fixture("seq.fa").to_str().unwrap(),
            "--sampler",
            "syncmer",
            "--zero",
        ])
        .run();

    // 4 sequences x 4 sequences = 16 pairs.
    assert_eq!(stdout.lines().count(), 16);
    // Self-comparison: identical syncmer set -> distance 0, jaccard 1.
    assert!(stdout.contains("seqA\tseqA\t0.0000\t1.0000\t1.0000"));
    // seqA (ACGT repeat) vs seqD (TGCA repeat): disjoint canonical syncmer sets.
    assert!(stdout.contains("seqA\tseqD\t1.0000\t0.0000\t0.0000"));
}

#[test]
fn command_dist_seq_syncmer_explicit_params() {
    // Explicit small s-mer/window so single-base differences fall in sampled
    // syncmers; seqA vs seqB differ only in the last base.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "seq",
            fixture("seq.fa").to_str().unwrap(),
            "--sampler",
            "syncmer",
            "-k",
            "5",
            "-w",
            "3",
            "--zero",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 16);
    // Self-comparison is always identity.
    assert!(stdout.contains("seqA\tseqA\t0.0000\t1.0000\t1.0000"));
}

#[test]
fn command_dist_seq_protein_rejects_mod_hasher() {
    // --hasher mod is DNA-only (canonical reverse complement); with --protein
    // it must error rather than silently producing garbage.
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "dist",
            "seq",
            fixture("seq.fa").to_str().unwrap(),
            "--hasher",
            "mod",
            "--protein",
        ])
        .run_fail();

    assert!(stderr.contains("--hasher mod is DNA-only"));
}

#[test]
fn command_dist_seq_syncmer_protein_defaults() {
    // --sampler syncmer --protein without explicit -k/-w applies protein
    // defaults (smer=7, window=5) rather than the DNA syng defaults (8, 55).
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "seq",
            fixture("seq.fa").to_str().unwrap(),
            "--sampler",
            "syncmer",
            "--protein",
            "--zero",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 16);
    // Self-comparison is identity regardless of alphabet interpretation.
    assert!(stdout.contains("seqA\tseqA\t0.0000\t1.0000\t1.0000"));
}
