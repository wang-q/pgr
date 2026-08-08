#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::path::PathBuf;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/dist/input")
        .join(name)
}

#[test]
fn command_dist_mash_basic() {
    // Mash-compatible bottom-k MinHash distances on random sequences.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "mash",
            fixture("random.fa").to_str().unwrap(),
            "-k",
            "21",
            "--size",
            "100",
            "--zero",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 9); // 3 x 3 pairs
    assert!(stdout.contains("r1\tr1\t0.0000\t1.0000\t1.0000"));
    assert!(stdout.contains("r1\tr2\t0.1221\t0.0400\t0.0800"));
    assert!(stdout.contains("r1\tr3\t1.0000\t0.0000\t0.0000"));
}

#[test]
fn command_dist_mash_merge() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "mash",
            fixture("random.fa").to_str().unwrap(),
            fixture("random.fa").to_str().unwrap(),
            "-k",
            "21",
            "--size",
            "100",
            "--merge",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 1);
    // Both files merge all sequences; identical merged sketches -> distance 0.
    assert!(stdout.contains("0.0000\t1.0000\t1.0000"));
}

#[test]
fn command_dist_mash_sim() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "mash",
            fixture("random.fa").to_str().unwrap(),
            "-k",
            "21",
            "--size",
            "100",
            "--zero",
            "--sim",
        ])
        .run();

    assert!(stdout.contains("r1\tr1\t1.0000\t1.0000\t1.0000"));
    // dist 0.1221 -> sim 0.8779
    assert!(stdout.contains("r1\tr2\t0.8779\t0.0400\t0.0800"));
}
