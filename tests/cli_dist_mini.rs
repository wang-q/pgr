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
fn command_dist_mini() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "mini",
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
fn command_dist_mini_sim() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "mini",
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
fn command_dist_mini_genome() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "mini",
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
fn command_dist_mini_merge() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "mini",
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
fn command_dist_mini_protein_rejects_mod_hasher() {
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "dist",
            "mini",
            fixture("seq.fa").to_str().unwrap(),
            "--hasher",
            "mod",
            "--protein",
        ])
        .run_fail();

    assert!(stderr.contains("--hasher mod is DNA-only"));
}
