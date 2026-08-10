#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::io::{Read, Write};

fn read_gz(path: &str) -> Vec<u8> {
    let mut out = Vec::new();
    let mut dec = flate2::read::GzDecoder::new(std::fs::File::open(path).unwrap());
    dec.read_to_end(&mut out).unwrap();
    out
}

fn write_temp(content: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::Builder::new().suffix(".fq").tempfile().unwrap();
    f.write_all(content.as_bytes()).unwrap();
    f
}

#[test]
fn command_fq_clump_matches_bbtools_clumpify_golden() {
    // Byte-level comparison against BBTools 39.38
    // `clumpify.sh seed=1` on the Lambda golden data (see
    // tests/bbtools/Lambda/README.md).
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("clump.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "clump",
            "tests/bbtools/Lambda/R1.fq.gz",
            "tests/bbtools/Lambda/R2.fq.gz",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read(&out).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/clumpify.fq.gz")
    );
}

#[test]
fn command_fq_clump_dedupe_removes_exact_pairs_only() {
    // Pair A appears twice (r1 low quality, r2 high quality); pair B shares
    // R1 with A but has a different R2 (not a duplicate); pair C is unique.
    // Dedupe must keep only the higher-quality copy of A, keep B, and keep C.
    let seq1 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let seq2 = "TGCAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
    let seq3 = "GGGGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
    let input = format!(
        "@r1/1 c1\n{seq1}\n+\n{}\n@r1/2 c2\n{seq2}\n+\n{}\n\
         @r2/1 c1\n{seq1}\n+\n{}\n@r2/2 c2\n{seq2}\n+\n{}\n\
         @r3/1 c1\n{seq1}\n+\n{}\n@r3/2 c2\n{seq3}\n+\n{}\n\
         @r4/1 c1\n{seq3}\n+\n{}\n@r4/2 c2\n{seq2}\n+\n{}\n",
        "!".repeat(seq1.len()),
        "!".repeat(seq2.len()),
        "I".repeat(seq1.len()),
        "I".repeat(seq2.len()),
        "I".repeat(seq1.len()),
        "I".repeat(seq3.len()),
        "I".repeat(seq3.len()),
        "I".repeat(seq2.len()),
    );
    let file = write_temp(&input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "clump",
            file.path().to_str().unwrap(),
            "--dedupe",
            "--dupesubs",
            "0",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    let out = std::fs::read_to_string(&out).unwrap();
    // Higher-quality copy of pair A kept; low-quality copy removed.
    assert!(out.contains("@r2/1"), "high-quality copy must survive");
    assert!(
        !out.contains("@r1/1"),
        "low-quality duplicate must be removed"
    );
    // B shares R1 with A but differs in R2: not a duplicate.
    assert!(out.contains("@r3/1"), "R1-only match must not be deduped");
    // Unique pair kept.
    assert!(out.contains("@r4/1"));
}

#[test]
fn command_fq_clump_external_bucket_path_matches_in_memory_set() {
    // --sort-mode bucket forces the external hash-bucket path. Order becomes
    // bucket-concatenated (documented), so compare read sets and determinism
    // rather than byte order; --buckets implies the same path with a fixed
    // count.
    let out_dir = tempfile::tempdir().unwrap();
    let mem = out_dir.path().join("mem.fq");
    let b1 = out_dir.path().join("bucket1.fq");
    let b2 = out_dir.path().join("bucket2.fq");
    let inputs = [
        "tests/bbtools/Lambda/R1.fq.gz",
        "tests/bbtools/Lambda/R2.fq.gz",
    ];

    PgrCmd::new()
        .args(&[
            "fq",
            "clump",
            inputs[0],
            inputs[1],
            "-o",
            mem.to_str().unwrap(),
        ])
        .assert()
        .success();
    for out in [&b1, &b2] {
        PgrCmd::new()
            .args(&[
                "fq",
                "clump",
                inputs[0],
                inputs[1],
                "--sort-mode",
                "bucket",
                "-o",
                out.to_str().unwrap(),
            ])
            .assert()
            .success();
    }
    PgrCmd::new()
        .args(&[
            "fq",
            "clump",
            inputs[0],
            inputs[1],
            "--buckets",
            "16",
            "-o",
            out_dir.path().join("b16.fq").to_str().unwrap(),
        ])
        .assert()
        .success();

    let read_set = |p: &std::path::Path| -> Vec<String> {
        std::fs::read_to_string(p)
            .unwrap()
            .lines()
            .step_by(4)
            .map(|l| l.split_whitespace().next().unwrap().to_string())
            .collect::<Vec<_>>()
    };
    // Deterministic bucket output.
    assert_eq!(read_set(&b1), read_set(&b2));
    // Same read set as the in-memory path (order differs).
    let mut mem_set = read_set(&mem);
    mem_set.sort();
    let mut b_set = read_set(&b1);
    b_set.sort();
    assert_eq!(mem_set, b_set);
    // Fixed bucket count is deterministic and set-equivalent too.
    assert_eq!(
        read_set(&b1).len(),
        read_set(&out_dir.path().join("b16.fq")).len()
    );
}

#[test]
fn command_fq_clump_groups_identical_pairs_together() {
    // Three pairs; two share an R1 sequence and sort by the same pivot
    // k-mer, so they stay adjacent.
    let mut input = String::new();
    for i in 0..3 {
        input.push_str(&format!(
            "@r{i}/1 c1\n{}\n+\n{}\n@r{i}/2 c2\n{}\n+\n{}\n",
            if i == 2 {
                "A".repeat(50)
            } else {
                "ACGTACGT".repeat(7)
            },
            "I".repeat(if i == 2 { 50 } else { 56 }),
            "G".repeat(50),
            "I".repeat(50),
        ));
    }
    let file = write_temp(&input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "clump",
            file.path().to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    let content = std::fs::read_to_string(&out).unwrap();
    let names: Vec<&str> = content
        .lines()
        .step_by(4)
        .map(|l| l.split_whitespace().next().unwrap())
        .collect();
    assert_eq!(names.len(), 6);
    let adjacent = names.windows(4).any(|w| {
        (w[0].starts_with("@r0") && w[2].starts_with("@r1"))
            || (w[0].starts_with("@r1") && w[2].starts_with("@r0"))
    });
    assert!(adjacent, "shared-kmer pairs must stay adjacent: {names:?}");
}
