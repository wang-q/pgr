#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::io::{Read, Write};

fn write_temp(content: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::Builder::new().suffix(".fq").tempfile().unwrap();
    f.write_all(content.as_bytes()).unwrap();
    f
}

fn read_gz(path: &str) -> Vec<u8> {
    let mut out = Vec::new();
    let mut dec = flate2::read::GzDecoder::new(std::fs::File::open(path).unwrap());
    dec.read_to_end(&mut out).unwrap();
    out
}

#[test]
fn command_fq_sample_matches_bbtools_reformat_golden() {
    // Byte-level comparison against BBTools 39.38
    // `reformat.sh samplebasestarget=1000000 sampleseed=1` on the Lambda
    // golden data (see tests/bbtools/Lambda/README.md).
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("sample.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "sample",
            "tests/bbtools/Lambda/golden/filter.fq.gz",
            "-o",
            out.to_str().unwrap(),
            "--bases",
            "1000000",
            "--seed",
            "1",
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read(&out).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/sample.fq.gz")
    );
}

#[test]
fn command_fq_sample_deterministic_and_targeted() {
    let mut input = String::new();
    for i in 0..100 {
        input.push_str(&format!(
            "@r{i}/1 c1\n{}\n+\n{}\n",
            "A".repeat(100),
            "I".repeat(100)
        ));
        input.push_str(&format!(
            "@r{i}/2 c2\n{}\n+\n{}\n",
            "T".repeat(100),
            "I".repeat(100)
        ));
    }
    let file = write_temp(&input);
    let out_dir = tempfile::tempdir().unwrap();
    let out1 = out_dir.path().join("s1.fq");
    let out2 = out_dir.path().join("s2.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "sample",
            file.path().to_str().unwrap(),
            "-o",
            out1.to_str().unwrap(),
            "--bases",
            "5000",
        ])
        .assert()
        .success();
    PgrCmd::new()
        .args(&[
            "fq",
            "sample",
            file.path().to_str().unwrap(),
            "-o",
            out2.to_str().unwrap(),
            "--bases",
            "5000",
        ])
        .assert()
        .success();

    let s1 = std::fs::read_to_string(&out1).unwrap();
    let s2 = std::fs::read_to_string(&out2).unwrap();
    assert_eq!(s1, s2, "same seed must give identical output");
    let n_reads = s1.lines().count() / 4;
    assert_eq!(n_reads, 50, "5000 bases at 100 bp per read");
    // Pairs are kept together: names alternate /1 /2.
    let names: Vec<&str> = s1
        .lines()
        .step_by(4)
        .map(|l| l.split_whitespace().next().unwrap())
        .collect();
    for pair in names.chunks(2) {
        assert!(pair[0].ends_with("/1") && pair[1].ends_with("/2"));
    }
}

#[test]
fn command_fq_sample_missing_bases_is_clap_error() {
    // Regression: omitting --bases used to panic on an unwrap of the missing
    // argument; it must now be a clean clap usage error (non-zero exit, no
    // panic).
    let file = write_temp("@r/1 c1\nAAAA\n+\nIIII\n");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "fq",
            "sample",
            file.path().to_str().unwrap(),
            "-o",
            "stdout",
        ])
        .run_fail();
    assert!(
        stderr.contains("--bases") || stderr.contains("required"),
        "stderr: {stderr}"
    );
}
