#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::io::Read;

fn read_gz(path: &str) -> Vec<u8> {
    let mut out = Vec::new();
    let mut dec = flate2::read::GzDecoder::new(std::fs::File::open(path).unwrap());
    dec.read_to_end(&mut out).unwrap();
    out
}

fn write_temp(content: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::Builder::new().suffix(".fq").tempfile().unwrap();
    std::io::Write::write_all(&mut f, content.as_bytes()).unwrap();
    f
}

#[test]
fn command_fq_filter_matches_bbtools_filter_golden() {
    // Byte-level comparison against BBTools 39.38
    // `bbduk.sh ... k=27 cardinality tossbrokenreads=t ordered=t` (filter
    // mode) on the Lambda golden data.
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("filter.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "filter",
            "tests/bbtools/Lambda/golden/trim.fq.gz",
            "--ref",
            "tests/bbtools/Lambda/illumina_adapters.fa",
            "--k",
            "27",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read(&out).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/filter.fq.gz")
    );
}

#[test]
fn command_fq_filter_stats_match_bbtools() {
    // Filter mode stats: no adapter kmers survive at k=27, so only headers.
    let out_dir = tempfile::tempdir().unwrap();
    let stats = out_dir.path().join("filter.stats.txt");

    PgrCmd::new()
        .args(&[
            "fq",
            "filter",
            "tests/bbtools/Lambda/golden/trim.fq.gz",
            "--ref",
            "tests/bbtools/Lambda/illumina_adapters.fa",
            "--k",
            "27",
            "--stats",
            stats.to_str().unwrap(),
            "-o",
            out_dir.path().join("out.fq").to_str().unwrap(),
        ])
        .assert()
        .success();

    let expected = concat!(
        "#File\ttests/bbtools/Lambda/golden/trim.fq.gz\n",
        "#Total\t36384\n",
        "#Matched\t0\t0.00000%\n",
        "#Name\tReads\tReadsPct\n",
    );
    assert_eq!(std::fs::read_to_string(&stats).unwrap(), expected);
}

#[test]
fn command_fq_filter_parallel_out_of_range_is_friendly_error() {
    // Regression: an out-of-range --parallel must be rejected with a friendly
    // error before a thread pool is created.
    let file = write_temp("@r1\nACGTACGT\n+\nIIIIIIII\n");
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "fq",
            "filter",
            file.path().to_str().unwrap(),
            "--ref",
            "tests/bbtools/Lambda/illumina_adapters.fa",
            "--parallel",
            "1000000",
            "-o",
            "stdout",
        ])
        .run_fail();
    assert!(
        stderr.contains("--parallel") || stderr.contains("1..=1024"),
        "stderr: {stderr}"
    );
}
