#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::io::Read;

fn read_gz(path: &str) -> Vec<u8> {
    let mut out = Vec::new();
    let mut dec = flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap());
    dec.read_to_end(&mut out).unwrap();
    out
}

/// Tadpole-compatible error correction + tossing over the first 2000 Lambda
/// pairs, byte-identical to BBTools 40.01 `assemble.Tadpole ... ecc tossjunk
/// tossdepth=2 tossuncorrectable threads=1` (tests/bbtools/Lambda/README.md).
#[test]
fn command_fq_ecc_matches_bbtools_golden() {
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("ecct.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "ecc",
            "tests/bbtools/Lambda/golden/ecco_sub.fq.gz",
            "-o",
            out.to_str().unwrap(),
            "--toss-junk",
            "--toss-depth",
            "2",
            "--toss-uncorrectable",
        ])
        .assert()
        .success();
    assert_eq!(
        std::fs::read(&out).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/ecct_sub.fq.gz")
    );
}

/// Without toss flags all reads are kept and only corrected.
#[test]
fn command_fq_ecc_keeps_all_without_toss() {
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("ecct.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "ecc",
            "tests/bbtools/Lambda/golden/ecco_sub.fq.gz",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();
    let out = std::fs::read(&out).unwrap();
    assert_eq!(std::str::from_utf8(&out).unwrap().lines().count(), 16000); // 4000 reads x 4 lines
}

/// Toss decisions are per-pair: a pair is dropped only when both mates fail.
#[test]
fn command_fq_ecc_tossdepth_drops_bad_pairs() {
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("ecct.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "ecc",
            "tests/bbtools/Lambda/golden/ecco_sub.fq.gz",
            "-o",
            out.to_str().unwrap(),
            "--toss-depth",
            "2",
        ])
        .assert()
        .success();
    let out = std::fs::read(&out).unwrap();
    assert_eq!(std::str::from_utf8(&out).unwrap().lines().count(), 11184); // matches BBTools tossdepth=2 subset
}

/// A zero k-mer length must fail cleanly instead of panicking.
#[test]
fn command_fq_ecc_rejects_zero_kmer() {
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "ecc",
            "tests/bbtools/Lambda/golden/ecco_sub.fq.gz",
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "0",
        ])
        .assert()
        .failure();
}

/// FASTA input has no quality scores; error detection uses the fixed quality
/// 20 (BBTools null-quality behavior) instead of indexing an empty vector.
#[test]
fn command_fq_ecc_handles_fasta_input() {
    let out_dir = tempfile::tempdir().unwrap();
    let infile = out_dir.path().join("in.fa");
    let out = out_dir.path().join("out.fa");
    std::fs::write(
        &infile,
        ">r1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n>r2\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC\n",
    )
    .unwrap();
    PgrCmd::new()
        .args(&[
            "fq",
            "ecc",
            infile.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();
    let out = std::fs::read_to_string(&out).unwrap();
    assert_eq!(out.lines().count(), 4);
}

/// `--parallel` values are validated but the pipeline stays single-pass.
#[test]
fn command_fq_ecc_parallel_is_validated() {
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "ecc",
            "tests/bbtools/Lambda/golden/ecco_sub.fq.gz",
            "-o",
            out.to_str().unwrap(),
            "--parallel",
            "abc",
        ])
        .assert()
        .failure();
    PgrCmd::new()
        .args(&[
            "fq",
            "ecc",
            "tests/bbtools/Lambda/golden/ecco_sub.fq.gz",
            "-o",
            out.to_str().unwrap(),
            "--parallel",
            "8",
        ])
        .assert()
        .success();
}
