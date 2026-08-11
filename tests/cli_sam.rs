#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

const SAM: &str = "\
@HD\tVN:1.6\tSO:unknown
@SQ\tSN:ut\tLN:1000
r1\t0\tut\t101\t255\t50M\t*\t0\t0\tACGT\tIIII
r2\t16\tut\t201\t255\t50M\t*\t0\t0\tACGT\tIIII
r3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII
";

/// Mapped records become 1-based inclusive ranges; unmapped records and
/// header lines are skipped.
#[test]
fn command_sam_to_rg_basic() {
    let out_dir = tempfile::tempdir().unwrap();
    let sam = out_dir.path().join("in.sam");
    std::fs::write(&sam, SAM).unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["sam", "to-rg", sam.to_str().unwrap()])
        .run();
    assert_eq!(stdout, "ut:101-150\nut:201-250\n");
}

/// The range spans the full reference-consuming CIGAR (M/D/N/=/X), with
/// non-consuming operations (I/S/H/P) excluded.
#[test]
fn command_sam_to_rg_cigar_span() {
    let out_dir = tempfile::tempdir().unwrap();
    let sam = out_dir.path().join("in.sam");
    std::fs::write(
        &sam,
        "r1\t0\tut\t11\t255\t10M2I8M3D15N5=2X\t*\t0\t0\tACGT\tIIII\n",
    )
    .unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["sam", "to-rg", sam.to_str().unwrap()])
        .run();
    // 10 + 8 + 3 + 15 + 5 + 2 = 43 reference bases.
    assert_eq!(stdout, "ut:11-53\n");
}

/// Reads SAM from stdin.
#[test]
fn command_sam_to_rg_stdin() {
    let (stdout, _) = PgrCmd::new()
        .args(&["sam", "to-rg", "stdin"])
        .stdin("r1\t0\tut\t5\t255\t10M\t*\t0\t0\tACGT\tIIII\n")
        .run();
    assert_eq!(stdout, "ut:5-14\n");
}

/// Malformed records are skipped by default and rejected with `--strict`.
#[test]
fn command_sam_to_rg_strict() {
    let out_dir = tempfile::tempdir().unwrap();
    let sam = out_dir.path().join("bad.sam");
    std::fs::write(
        &sam,
        "r1\t0\tut\t5\t255\t10M\t*\t0\t0\tACGT\tIIII\nonly_two_fields\n",
    )
    .unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["sam", "to-rg", sam.to_str().unwrap()])
        .run();
    assert_eq!(stdout, "ut:5-14\n");
    let (_, stderr) = PgrCmd::new()
        .args(&["sam", "to-rg", sam.to_str().unwrap(), "--strict"])
        .run_fail();
    assert!(stderr.contains("FLAG"), "stderr: {stderr}");
}
