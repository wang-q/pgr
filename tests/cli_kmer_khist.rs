#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

#[test]
fn command_kmer_hist_kmercountexact_text_matches_bbtools_golden() {
    // Byte-level comparison against BBTools 39.38
    // `kmercountexact.sh khist=R.khist.txt peaks=R.peaks.txt k=31` on the
    // Lambda golden data (see tests/bbtools/Lambda/README.md).
    let out_dir = tempfile::tempdir().unwrap();
    let hist = out_dir.path().join("h.hist");
    let khist = out_dir.path().join("khist.txt");
    let peaks = out_dir.path().join("peaks.txt");

    PgrCmd::new()
        .args(&[
            "kmer",
            "hist",
            "tests/bbtools/Lambda/golden/filter.fq.gz",
            "-k",
            "31",
            "-o",
            hist.to_str().unwrap(),
            "--khist-text",
            khist.to_str().unwrap(),
            "--peaks",
            peaks.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read(&khist).unwrap(),
        std::fs::read("tests/bbtools/Lambda/golden/R.khist.txt").unwrap()
    );
    assert_eq!(
        std::fs::read(&peaks).unwrap(),
        std::fs::read("tests/bbtools/Lambda/golden/R.peaks.txt").unwrap()
    );
}

#[test]
fn command_kmer_hist_kmercountexact_works_from_pkt_table() {
    let out_dir = tempfile::tempdir().unwrap();
    let pkt = out_dir.path().join("t.pkt");
    let hist = out_dir.path().join("h.hist");
    let peaks = out_dir.path().join("peaks.txt");

    PgrCmd::new()
        .args(&[
            "kmer",
            "table",
            "tests/bbtools/Lambda/golden/filter.fq.gz",
            "-k",
            "31",
            "-o",
            pkt.to_str().unwrap(),
        ])
        .assert()
        .success();
    PgrCmd::new()
        .args(&[
            "kmer",
            "hist",
            "--table",
            pkt.to_str().unwrap(),
            "-o",
            hist.to_str().unwrap(),
            "--peaks",
            peaks.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read(&peaks).unwrap(),
        std::fs::read("tests/bbtools/Lambda/golden/R.peaks.txt").unwrap()
    );
}
