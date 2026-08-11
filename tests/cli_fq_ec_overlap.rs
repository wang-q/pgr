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

fn lambda(args: &[&str], out: &str, outu: Option<&str>, ihist: Option<&str>) {
    let mut full: Vec<&str> = vec![
        "fq",
        "ec-overlap",
        "tests/bbtools/Lambda/R1.2k.fq.gz",
        "tests/bbtools/Lambda/R2.2k.fq.gz",
        "-o",
        out,
    ];
    if let Some(u) = outu {
        full.push("--outu");
        full.push(u);
    }
    if let Some(h) = ihist {
        full.push("--ihist");
        full.push(h);
    }
    full.extend_from_slice(args);
    PgrCmd::new().args(&full).assert().success();
}

#[test]
fn command_fq_ec_overlap_matches_bbtools_golden() {
    // BBTools 40.01 `bbmerge.sh ... ecco mix vstrict` with the bundled
    // bbmerge.bbnet overlap filter, ordered, threads=1 (see
    // tests/bbtools/Lambda/README.md merge section).
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("ecco.fq");
    let ihist = out_dir.path().join("ihist1.txt");
    lambda(
        &[
            "--vstrict",
            "--net",
            "tests/bbtools/Lambda/golden/bbmerge.bbnet",
        ],
        out.to_str().unwrap(),
        None,
        Some(ihist.to_str().unwrap()),
    );
    assert_eq!(
        std::fs::read(&out).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/merge.ecco.fq.gz")
    );
    assert_eq!(
        std::fs::read(&ihist).unwrap(),
        std::fs::read("tests/bbtools/Lambda/golden/merge.ihist1.txt").unwrap()
    );
}

#[test]
fn command_fq_ec_overlap_novector_matches_bbtools_golden() {
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("ecco.fq");
    lambda(
        &["--vstrict", "--no-make-vector"],
        out.to_str().unwrap(),
        None,
        None,
    );
    assert_eq!(
        std::fs::read(&out).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/merge.novector.ecco.fq.gz")
    );
}

#[test]
fn command_fq_ec_overlap_without_mix_keeps_all_reads() {
    // `bbmerge.sh ... ecco` without an explicit `mix` auto-sets
    // MIX_BAD_AND_GOOD, so every pair is written to the main output; the
    // result must match the golden that was generated with `ecco mix`.
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("ecco.fq");
    lambda(
        &[
            "--vstrict",
            "--net",
            "tests/bbtools/Lambda/golden/bbmerge.bbnet",
        ],
        out.to_str().unwrap(),
        None,
        None,
    );
    assert_eq!(
        std::fs::read(&out).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/merge.ecco.fq.gz")
    );
}

#[test]
fn command_fq_ec_overlap_no_mix_writes_only_corrected_pairs() {
    // `bbmerge ... ecco mix=f`: only overlapping pairs are corrected; the
    // rest are dropped when no --outu is given.
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("ecco.fq");
    lambda(
        &[
            "--no-mix",
            "--vstrict",
            "--net",
            "tests/bbtools/Lambda/golden/bbmerge.bbnet",
        ],
        out.to_str().unwrap(),
        None,
        None,
    );
    let out = std::fs::read(&out).unwrap();
    assert_eq!(std::str::from_utf8(&out).unwrap().lines().count(), 80);
}
