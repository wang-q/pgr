#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;

/// Deterministic pseudo-random DNA (LCG, no ACGT periodicity).
fn random_seq(len: usize, seed: u64) -> String {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut x = seed;
    (0..len)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            bases[(x >> 33) as usize & 3] as char
        })
        .collect()
}

/// Omitting the query implies self-alignment: with a multi-file target
/// directory, only same-basename pairs are aligned (not the full n x n
/// cross product). Regression for the `--self` / omitted-query mismatch.
#[test]
fn command_align_lastz_omitted_query_is_self() {
    if which::which("lastz").is_err() {
        eprintln!("Skipping command_align_lastz_omitted_query_is_self: lastz not installed");
        return;
    }
    let temp = tempfile::TempDir::new().unwrap();
    let target = temp.path().join("targets");
    fs::create_dir_all(&target).unwrap();
    fs::write(
        target.join("a.fa"),
        format!(">a\n{}\n{}\n", random_seq(1200, 11), random_seq(400, 12)),
    )
    .unwrap();
    fs::write(
        target.join("b.fa"),
        format!(">b\n{}\n{}\n", random_seq(1200, 13), random_seq(400, 14)),
    )
    .unwrap();

    let out = temp.path().join("out");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "lastz",
            target.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(
        !stderr.contains("error:"),
        "self alignment failed: {stderr}"
    );

    let lavs: Vec<String> = fs::read_dir(&out)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().to_string())
        .filter(|n| n.ends_with(".lav"))
        .collect();
    assert_eq!(lavs.len(), 2, "expected only self pairs, got: {lavs:?}");
    assert!(
        lavs.iter().any(|n| n.contains("[a]vs[a]")),
        "missing [a]vs[a].lav: {lavs:?}"
    );
    assert!(
        lavs.iter().any(|n| n.contains("[b]vs[b]")),
        "missing [b]vs[b].lav: {lavs:?}"
    );
}

/// Self mode on a directory whose files share a basename must align each file
/// to itself only; two files that merely share a name (e.g. `a/dup.fa` and
/// `b/dup.fa`) must not be cross-aligned (the old basename comparison let
/// them through as spurious pairs).
#[test]
fn command_align_lastz_self_duplicate_basenames() {
    if which::which("lastz").is_err() {
        eprintln!("Skipping command_align_lastz_self_duplicate_basenames: lastz not installed");
        return;
    }
    let temp = tempfile::TempDir::new().unwrap();
    let a = temp.path().join("a");
    let b = temp.path().join("b");
    fs::create_dir_all(&a).unwrap();
    fs::create_dir_all(&b).unwrap();
    fs::write(a.join("dup.fa"), format!(">a\n{}\n", random_seq(1200, 21))).unwrap();
    fs::write(b.join("dup.fa"), format!(">b\n{}\n", random_seq(1200, 22))).unwrap();

    let out = temp.path().join("out");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "align",
            "lastz",
            temp.path().to_str().unwrap(),
            "--self",
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert!(
        !stderr.contains("error:"),
        "self alignment failed: {stderr}"
    );

    let lavs: Vec<String> = fs::read_dir(&out)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().to_string())
        .filter(|n| n.ends_with(".lav"))
        .collect();
    assert_eq!(
        lavs.len(),
        2,
        "expected exactly the two self pairs, got: {lavs:?}"
    );
}
