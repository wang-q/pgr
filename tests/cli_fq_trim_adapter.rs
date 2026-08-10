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
fn command_fq_trim_adapter_matches_bbtools_trim_golden() {
    // Byte-level comparison against BBTools 39.38
    // `bbduk.sh ... ktrim=r k=23 mink=11 hdist=1 tbo tpe qtrim=r trimq=15
    // minlen=60 maxns=0 ftm=5 tossbrokenreads=t ordered=t` on the Lambda
    // golden data (see tests/bbtools/Lambda/README.md).
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("trim.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "trim-adapter",
            "tests/bbtools/Lambda/golden/clumpify.fq.gz",
            "--ref",
            "tests/bbtools/Lambda/illumina_adapters.fa",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read(&out).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/trim.fq.gz")
    );
}

#[test]
fn command_fq_trim_adapter_filter_matches_bbtools_filter_golden() {
    // Byte-level comparison against BBTools 39.38
    // `bbduk.sh ... k=27 cardinality tossbrokenreads=t ordered=t` (filter
    // mode) on the Lambda golden data.
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("filter.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "trim-adapter",
            "tests/bbtools/Lambda/golden/trim.fq.gz",
            "--ref",
            "tests/bbtools/Lambda/illumina_adapters.fa",
            "--no-ktrim",
            "--no-tbo",
            "--no-tpe",
            "--no-qtrim",
            "--k",
            "27",
            "--mink",
            "0",
            "--minlen",
            "0",
            "--maxns=-1",
            "--ftm",
            "0",
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
fn command_fq_trim_adapter_removes_adapter_and_keeps_clean_read() {
    // A read whose 3' end is a known adapter is trimmed; a clean read is
    // untouched except the ftm multiple-of-5 normalization.
    let input = format!(
        "@r1\n{}AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n+\n{}IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        "A".repeat(70),
        "I".repeat(70)
    );
    let file = write_temp(&input);
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "trim-adapter",
            file.path().to_str().unwrap(),
            "--ref",
            "tests/bbtools/Lambda/illumina_adapters.fa",
            "--no-tbo",
            "--no-tpe",
            "--maxns=-1",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    let out = std::fs::read_to_string(&out).unwrap();
    let seq = out.lines().nth(1).unwrap();
    assert!(
        !seq.contains("AATGATACGGCGACCACC"),
        "adapter must be trimmed"
    );
    // The hdist=1 table can cut one extra base at the adapter boundary, but
    // the read must keep well above the 60 bp minimum.
    assert!(seq.len() >= 60, "prefix must survive: {seq}");
    assert_eq!(seq.len(), 69, "bbduk-compatible cut position");
}
