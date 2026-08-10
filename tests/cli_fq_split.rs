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
fn command_fq_split_interleaved_into_r1_r2_and_singles() {
    let input = "\
@r1/1 c1
ACGT
+
!!!!
@r1/2 c2
TGCA
+
####
@solo/1 c3
AAAA
+
BBBB
";
    let file = write_temp(input);
    let out_dir = tempfile::tempdir().unwrap();
    let r1 = out_dir.path().join("r1.fq");
    let r2 = out_dir.path().join("r2.fq");
    let s = out_dir.path().join("s.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "split",
            file.path().to_str().unwrap(),
            "-o",
            r1.to_str().unwrap(),
            "--outfile-2",
            r2.to_str().unwrap(),
            "--outfile-single",
            s.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert_eq!(
        std::fs::read_to_string(&r1).unwrap(),
        "@r1/1 c1\nACGT\n+\n!!!!\n"
    );
    assert_eq!(
        std::fs::read_to_string(&r2).unwrap(),
        "@r1/2 c2\nTGCA\n+\n####\n"
    );
    assert_eq!(
        std::fs::read_to_string(&s).unwrap(),
        "@solo/1 c3\nAAAA\n+\nBBBB\n"
    );
}

#[test]
fn command_fq_split_matches_bbtools_repair_golden() {
    // Byte-level comparison against BBTools 39.38 `repair.sh rp` output on the
    // Lambda golden data (see tests/bbtools/Lambda/README.md).
    let out_dir = tempfile::tempdir().unwrap();
    let r1 = out_dir.path().join("r1.fq");
    let r2 = out_dir.path().join("r2.fq");
    let s = out_dir.path().join("s.fq");

    PgrCmd::new()
        .args(&[
            "fq",
            "split",
            "tests/bbtools/Lambda/golden/filter.fq.gz",
            "-o",
            r1.to_str().unwrap(),
            "--outfile-2",
            r2.to_str().unwrap(),
            "--outfile-single",
            s.to_str().unwrap(),
        ])
        .assert()
        .success();

    for (name, golden) in [
        ("r1", "tests/bbtools/Lambda/golden/R1.fq.gz"),
        ("r2", "tests/bbtools/Lambda/golden/R2.fq.gz"),
        ("s", "tests/bbtools/Lambda/golden/Rs.fq.gz"),
    ] {
        let path = match name {
            "r1" => &r1,
            "r2" => &r2,
            _ => &s,
        };
        assert_eq!(
            std::fs::read(path).unwrap(),
            read_gz(golden),
            "{name} differs from golden"
        );
    }
}
