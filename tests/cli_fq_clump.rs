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
