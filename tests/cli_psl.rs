#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

fn get_path(subcommand: &str, dir: &str, filename: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("tests/psl");
    path.push(subcommand);
    path.push(dir);
    path.push(filename);
    path
}

//
// psl histo
//

#[test]
fn test_histo_apq_base() {
    let temp = TempDir::new().unwrap();
    let input = get_path("histo", "input", "basic.psl");
    let output = temp.path().join("apq.histo");

    PgrCmd::new()
        .args(&[
            "psl",
            "histo",
            "--what",
            "alignsPerQuery",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    // Check for expected counts.
    // NM_033178.1: 2
    // NM_173571.1: 3
    // NM_000014.3: 1
    // NM_000015.1: 1
    // NM_153248.2: 1
    // NM_005577.1: 4
    // NM_FAKE.1: 2
    // Expected output order depends on hash map iteration unless sorted.
    // I implemented sorting by key.
    // Sorted keys: NM_000014.3, NM_000015.1, NM_005577.1, NM_033178.1, NM_153248.2, NM_173571.1, NM_FAKE.1
    // Counts: 1, 1, 4, 2, 1, 3, 2
    let expected = "1\n1\n4\n2\n1\n3\n2\n";
    assert_eq!(output_content, expected);
}

#[test]
fn test_histo_apq_multi() {
    let temp = TempDir::new().unwrap();
    let input = get_path("histo", "input", "basic.psl");
    let output = temp.path().join("apq_multi.histo");

    PgrCmd::new()
        .args(&[
            "psl",
            "histo",
            "--what",
            "alignsPerQuery",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--multi-only",
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    // Multi only: NM_005577.1 (4), NM_033178.1 (2), NM_173571.1 (3), NM_FAKE.1 (2)
    // Order: NM_005577.1, NM_033178.1, NM_173571.1, NM_FAKE.1
    // Counts: 4, 2, 3, 2
    let expected = "4\n2\n3\n2\n";
    assert_eq!(output_content, expected);
}

#[test]
fn test_histo_cover_spread() {
    let temp = TempDir::new().unwrap();
    let input = get_path("histo", "input", "basic.psl");
    let output = temp.path().join("cover.histo");

    PgrCmd::new()
        .args(&[
            "psl",
            "histo",
            "--what",
            "coverSpread",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    // NM_000014.3: 1 align. Spread = 0.
    // NM_005577.1: 4 aligns.
    //   3335+96+0 / 13938 = 0.24616
    //   3444+105+0 / 13938 = 0.25463
    //   3482+120+0 / 13938 = 0.25843
    //   6410+4+0 / 13938 = 0.46018
    //   Diff: 0.46018 - 0.24616 = 0.2140
    //
    // Just checking it runs and produces output. Precise float matching is tricky.
    // I will check if output contains "0.2140"
    let output_content = fs::read_to_string(&output).unwrap();
    assert!(output_content.contains("0.2140"));
    assert!(output_content.contains("0.0000")); // Singletons or identicals
}

#[test]
fn test_histo_id_spread() {
    let temp = TempDir::new().unwrap();
    let input = get_path("histo", "input", "basic.psl");
    let output = temp.path().join("id.histo");

    PgrCmd::new()
        .args(&[
            "psl",
            "histo",
            "--what",
            "idSpread",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    let lines: Vec<&str> = output_content.lines().collect();
    // basic.psl has 7 unique queries.
    assert_eq!(lines.len(), 7);
}

//
// psl to-chain
//

#[test]
fn test_to_chain_fix_strand() {
    let temp = TempDir::new().unwrap();
    let input = get_path("to_chain", "input", "mtor.psl");
    let expected_output = get_path("to_chain", "expected", "example3.chain");
    let output = temp.path().join("out.chain");

    PgrCmd::new()
        .args(&[
            "psl",
            "to-chain",
            input.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--fix-strand",
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    let expected_content = fs::read_to_string(&expected_output).unwrap();
    assert_eq!(output_content, expected_content);
}

#[test]
fn test_to_chain_fail_neg_strand() {
    let temp = TempDir::new().unwrap();
    let input = get_path("to_chain", "input", "mtor.psl");
    let output = temp.path().join("out.chain");

    PgrCmd::new()
        .args(&[
            "psl",
            "to-chain",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run_fail();
}

//
// psl rc
//

#[test]
fn test_rc_mrna() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "psl",
            "rc",
            get_path("rc", "input", "mrna.psl").to_str().unwrap(),
            "-o",
            "stdout",
        ])
        .run();

    let expected = std::fs::read_to_string(get_path("rc", "expected", "mrnaTest.psl")).unwrap();
    assert_eq!(stdout.replace("\r\n", "\n"), expected.replace("\r\n", "\n"));
}

#[test]
fn test_rc_trans() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "psl",
            "rc",
            get_path("rc", "input", "trans.psl").to_str().unwrap(),
            "-o",
            "stdout",
        ])
        .run();

    let expected = std::fs::read_to_string(get_path("rc", "expected", "transTest.psl")).unwrap();
    assert_eq!(stdout.replace("\r\n", "\n"), expected.replace("\r\n", "\n"));
}

//
// psl lift
//

#[test]
fn test_lift_basic() {
    let temp = TempDir::new().unwrap();
    let input = get_path("lift", "", "test_fragment.psl");
    let sizes = get_path("lift", "", "chrom.sizes");
    let output = temp.path().join("lifted.psl");

    PgrCmd::new()
        .args(&[
            "psl",
            "lift",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--q-sizes",
            sizes.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();

    // Expected output check
    // The input file contains two records.
    // First record: chr1:101-200 (+), qStart=10, qEnd=20 (on fragment).
    //   Lifted: chr1 (+), qStart=100+10=110, qEnd=100+20=120.
    // Second record: chr1:101-200 (-), qStart=10, qEnd=20 (on RC fragment).
    //   Lifted: chr1 (-), qStart=810, qEnd=820.
    //   Calculation for (-):
    //     q_size=1000. offset=100.
    //     new_qStart = 1000 - (100 - 10 + 100) = 810.
    //     new_qEnd = 1000 - (100 - 20 + 100) = 820.

    // Check first record
    assert!(output_content.contains("chr1\t1000\t110\t120"));
    // Check second record
    assert!(output_content.contains("chr1\t1000\t810\t820"));
    // Check that qStarts for blocks are also correct
    // First record block start: 110
    assert!(output_content.contains("110,\t500,"));
    // Second record block start: 810
    assert!(output_content.contains("810,\t500,"));
}

#[test]
fn test_lift_target() {
    let temp = TempDir::new().unwrap();
    let input = get_path("lift", "", "target_lift.psl");
    let sizes = get_path("lift", "", "chrom.sizes");
    let output = temp.path().join("target_lifted.psl");

    PgrCmd::new()
        .args(&[
            "psl",
            "lift",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--t-sizes",
            sizes.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();

    // Expected output check
    // First record: Target chr1:101-200 (+). tStart=10, tEnd=60.
    //   Lifted: chr1 (+). tStart=110, tEnd=160.
    // Second record: Target chr1:101-200 (-). tStart=10, tEnd=60.
    //   Lifted: chr1 (-). tStart=810, tEnd=860. (Size 1000)

    // Check first record
    // Target name: chr1
    // Target size: 1000
    // Target start: 110
    // Target end: 160
    assert!(output_content.contains("seq1\t100\t0\t50\tchr1\t1000\t110\t160"));

    // Check second record
    // Target name: chr1
    // Target size: 1000
    // Target start: 810
    // Target end: 860
    assert!(output_content.contains("seq1\t100\t0\t50\tchr1\t1000\t810\t860"));

    // Check tStarts
    // First: 110
    assert!(output_content.contains(",\t110,"));
    // Second: 810
    assert!(output_content.contains(",\t810,"));
}

#[test]
fn test_lift_fail() {
    // Missing arguments
    PgrCmd::new().args(&["psl", "lift"]).run_fail();
}

//
// psl stats
//

#[test]
fn test_stats_basic() {
    let temp = TempDir::new().unwrap();
    let input = get_path("stats", "input", "stats_basic.psl");
    let output = temp.path().join("stats.tsv");

    PgrCmd::new()
        .args(&[
            "psl",
            "stats",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    let lines: Vec<&str> = output_content.lines().collect();
    // Default is per-alignment stats.
    // Input has 31 records. Output should have 32 lines (header + 31).
    assert_eq!(lines.len(), 32);
}

//
// psl to-range
//

#[test]
fn test_to_range_basic() {
    let temp = TempDir::new().unwrap();
    let input = get_path("lift", "", "test_fragment.psl");
    let output = temp.path().join("ranges.rg");

    PgrCmd::new()
        .args(&[
            "psl",
            "to-range",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();

    // Check output content
    // Input:
    // 1. chr1:101-200 (+), qStart=10, qEnd=20.
    //    Range: chr1:101-200:11-20
    // 2. chr1:101-200 (-), qStart=10, qEnd=20.
    //    qSize=100.
    //    Range: chr1:101-200:81-90 (as calculated before)

    let lines: Vec<&str> = output_content.lines().collect();
    assert_eq!(lines.len(), 2);
    assert_eq!(lines[0], "chr1:101-200:11-20");
    assert_eq!(lines[1], "chr1:101-200:81-90");
}
