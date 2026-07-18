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
