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

#[test]
fn command_fq_split_stdout_matches_golden() {
    // R1 written to stdout must match the BBTools repair golden byte for
    // byte; R2/singles files match as in the file-output test.
    let out_dir = tempfile::tempdir().unwrap();
    let r2 = out_dir.path().join("r2.fq");
    let s = out_dir.path().join("s.fq");

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fq",
            "split",
            "tests/bbtools/Lambda/golden/filter.fq.gz",
            "-o",
            "stdout",
            "--outfile-2",
            r2.to_str().unwrap(),
            "--outfile-single",
            s.to_str().unwrap(),
        ])
        .run();

    assert_eq!(
        stdout.as_bytes(),
        read_gz("tests/bbtools/Lambda/golden/R1.fq.gz"),
        "stdout R1 differs from golden"
    );
    assert_eq!(
        std::fs::read(&r2).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/R2.fq.gz"),
        "R2 differs from golden"
    );
    assert_eq!(
        std::fs::read(&s).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/Rs.fq.gz"),
        "singles differs from golden"
    );
}

#[test]
fn command_fq_split_without_singles_discards_trailing_read() {
    // A trailing read without its mate is discarded with a warning when no
    // --outfile-single is given.
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

    let (_, stderr) = PgrCmd::new()
        .args(&[
            "fq",
            "split",
            file.path().to_str().unwrap(),
            "-o",
            r1.to_str().unwrap(),
            "--outfile-2",
            r2.to_str().unwrap(),
        ])
        .run();

    assert!(
        stderr.contains("unpaired read discarded"),
        "stderr: {stderr}"
    );
    assert_eq!(
        std::fs::read_to_string(&r1).unwrap(),
        "@r1/1 c1\nACGT\n+\n!!!!\n"
    );
    assert_eq!(
        std::fs::read_to_string(&r2).unwrap(),
        "@r1/2 c2\nTGCA\n+\n####\n"
    );
}

#[test]
fn command_fq_split_repair_pairs_by_name() {
    // repair.sh `rp` mode: r2/2 comes before r2/1 (disordered), r3/1 is an
    // orphan, and "orphan" has no pair marker. --repair must recover the
    // order and route singletons to --outfile-single.
    let input = "\
@r1/1
AAAA
+
IIII
@r1/2
TTTT
+
IIII
@r2/2
CCCC
+
IIII
@r2/1
GGGG
+
IIII
@r3/1
ACAC
+
IIII
@orphan
GTGT
+
IIII
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
            "--repair",
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
        "@r1/1\nAAAA\n+\nIIII\n@r2/1\nGGGG\n+\nIIII\n"
    );
    assert_eq!(
        std::fs::read_to_string(&r2).unwrap(),
        "@r1/2\nTTTT\n+\nIIII\n@r2/2\nCCCC\n+\nIIII\n"
    );
    assert_eq!(
        std::fs::read_to_string(&s).unwrap(),
        "@r3/1\nACAC\n+\nIIII\n@orphan\nGTGT\n+\nIIII\n"
    );
}

#[test]
fn command_fq_split_repair_matches_bbtools_golden() {
    // On well-ordered input (BBTools filter golden), --repair must produce
    // the same R1/R2/Rs as repair.sh rp mode.
    let out_dir = tempfile::tempdir().unwrap();
    let r1 = out_dir.path().join("r1.fq");
    let r2 = out_dir.path().join("r2.fq");
    let s = out_dir.path().join("s.fq");
    PgrCmd::new()
        .args(&[
            "fq",
            "split",
            "tests/bbtools/Lambda/golden/filter.fq.gz",
            "--repair",
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
        std::fs::read(&r1).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/R1.fq.gz")
    );
    assert_eq!(
        std::fs::read(&r2).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/R2.fq.gz")
    );
    assert_eq!(
        std::fs::read(&s).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/Rs.fq.gz")
    );
}

#[test]
fn command_fq_split_missing_outfile2_is_clap_error() {
    // Regression: omitting --outfile-2 used to panic on an unwrap of the
    // missing argument; it must now be a clean clap usage error (non-zero
    // exit, no panic).
    let input = "\
@r1/1 c1
ACGT
+
!!!!
@r1/2 c2
TGCA
+
####
";
    let file = write_temp(input);
    let (_, stderr) = PgrCmd::new()
        .args(&["fq", "split", file.path().to_str().unwrap(), "-o", "stdout"])
        .run_fail();
    assert!(
        stderr.contains("--outfile-2") || stderr.contains("required"),
        "stderr: {stderr}"
    );
}
