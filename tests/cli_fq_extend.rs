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

/// Tadpole-compatible read extension (k=62, el=20, er=20) over the first
/// 2000 Lambda pairs, byte-identical to BBTools 40.01 `assemble.Tadpole ...
/// mode=extend el=20 er=20 k=62 threads=1` (Tadpole2 long-k-mer path).
#[test]
fn command_fq_extend_matches_bbtools_golden() {
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("ext.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "extend",
            "tests/bbtools/Lambda/golden/ecco_sub.fq.gz",
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "62",
            "--el",
            "20",
            "--er",
            "20",
        ])
        .assert()
        .success();
    assert_eq!(
        std::fs::read(&out).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/ext_sub.fq.gz")
    );
}

/// Defaults (k=31) also run; reads are never discarded in extend mode.
#[test]
fn command_fq_extend_defaults_keep_all_reads() {
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("ext.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "extend",
            "tests/bbtools/Lambda/golden/ecco_sub.fq.gz",
            "-o",
            out.to_str().unwrap(),
            "--el",
            "10",
            "--er",
            "10",
        ])
        .assert()
        .success();
    let out = std::fs::read(&out).unwrap();
    assert_eq!(std::str::from_utf8(&out).unwrap().lines().count(), 16000);
}

/// A zero k-mer length must fail cleanly instead of panicking.
#[test]
fn command_fq_extend_rejects_zero_kmer() {
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "extend",
            "tests/bbtools/Lambda/golden/ecco_sub.fq.gz",
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "0",
        ])
        .assert()
        .failure();
}
