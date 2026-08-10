#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::io::Write;

fn write_temp(content: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::Builder::new().suffix(".fq").tempfile().unwrap();
    f.write_all(content.as_bytes()).unwrap();
    f
}

#[test]
fn command_fq_norm_keeps_high_depth_reads_and_tosses_low_depth() {
    // A 100x-covered read survives; a read with unique (depth-1) k-mers is
    // tossed.
    let mut input = String::new();
    for i in 0..100 {
        input.push_str(&format!(
            "@hi{i}\n{}\n+\n{}\n",
            "ACGT".repeat(20),
            "I".repeat(80)
        ));
    }
    // A random 64-mer whose 31-mers are unique in the dataset (depth 1).
    input.push_str(&format!(
        "@lo\n{}\n+\n{}\n",
        "GATCCTAGACGTTCGATCGGTACCTAGCATGCAGTTACGTACGATCGTAGCTAGCGGATCGATC",
        "I".repeat(64)
    ));
    let file = write_temp(&input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "norm",
            file.path().to_str().unwrap(),
            "-k",
            "31",
            "--min",
            "3",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    let out = std::fs::read_to_string(&out).unwrap();
    assert!(out.contains("@hi99"));
    assert!(!out.contains("@lo"));
}
