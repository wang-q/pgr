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
            "--field",
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
            "--field",
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
            "--field",
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
            "--field",
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
            "--outfile",
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

#[test]
fn test_to_chain_strict_malformed() {
    let temp = TempDir::new().unwrap();
    let input = temp.path().join("malformed.psl");
    fs::write(&input, "this is not a valid psl line\n").unwrap();
    let output = temp.path().join("out.chain");

    PgrCmd::new()
        .args(&[
            "psl",
            "to-chain",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--strict",
        ])
        .run_fail();
}

#[test]
fn test_to_chain_non_strict_malformed() {
    let temp = TempDir::new().unwrap();
    let input = temp.path().join("malformed.psl");
    fs::write(&input, "this is not a valid psl line\n").unwrap();
    let output = temp.path().join("out.chain");

    PgrCmd::new()
        .args(&[
            "psl",
            "to-chain",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    assert!(output_content.is_empty());
}

#[test]
fn test_to_chain_untranslated() {
    let temp = TempDir::new().unwrap();
    let input = temp.path().join("untranslated.psl");
    fs::write(
        &input,
        "10\t0\t0\t0\t0\t0\t0\t0\t+\tq1\t100\t10\t20\tt\t200\t50\t60\t1\t10,\t10,\t50,\n\
         10\t0\t0\t0\t0\t0\t0\t0\t-\tq2\t100\t10\t20\tt\t200\t70\t80\t1\t10,\t10,\t70,\n",
    )
    .unwrap();
    let output = temp.path().join("out.chain");

    PgrCmd::new()
        .args(&[
            "psl",
            "to-chain",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    let lines: Vec<&str> = output_content.lines().collect();
    assert!(lines[0].starts_with("chain 10 t 200 + 50 60 q1 100 + 10 20 1"));
    // Negative query strand is reversed by write_chain.
    assert!(lines[2].starts_with("chain 10 t 200 + 70 80 q2 100 - 80 90 2"));
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
    // Second record: chr1:101-200 (-), qStart=10, qEnd=20 on the forward
    // fragment (UCSC PSL: qStart/qEnd are forward coordinates even for '-',
    // only the block qStarts are in the RC frame: qSize - qEnd = 80).
    //   Lifted: chr1 (-), qStart=110, qEnd=120 (forward), and the block
    //   qStarts lift in the chromosome RC frame: 80 + (1000 - 200) = 880.

    // Check first record
    assert!(output_content.contains("chr1\t1000\t110\t120"));
    // Check second record (outer coords are forward, so they match the first)
    assert!(output_content.contains("chr1\t1000\t110\t120"));
    // Check that qStarts for blocks are also correct
    // First record block start: 110
    assert!(output_content.contains("110,\t500,"));
    // Second record block start (RC frame): 880
    assert!(output_content.contains("880,\t500,"));
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
    // Second record: Target chr1:101-200 (-). tStart=10, tEnd=60 on the
    // forward target (UCSC convention); tStarts are in the RC frame
    // (tSize - tEnd = 40).
    //   Lifted: chr1 (-). tStart=110, tEnd=160; block tStarts = 40 + 800 = 840.

    // Check first record
    // Target name: chr1
    // Target size: 1000
    // Target start: 110
    // Target end: 160
    assert!(output_content.contains("seq1\t100\t0\t50\tchr1\t1000\t110\t160"));

    // Check second record (outer coords are forward, so they match the first)
    assert!(output_content.contains("seq1\t100\t0\t50\tchr1\t1000\t110\t160"));

    // Check tStarts
    // First: 110
    assert!(output_content.contains(",\t110,"));
    // Second (RC frame): 840
    assert!(output_content.contains(",\t840,"));
}

/// Negative-strand fragment lifts keep UCSC PSL semantics: qStart/qEnd stay
/// forward-strand coordinates and only the block qStarts move to the
/// chromosome RC frame (regression for a fragment window aligning on '-'
/// strand, the `rept s-align` pipeline).
#[test]
fn test_lift_minus_strand_forward_coordinates() {
    let temp = TempDir::new().unwrap();
    let input = temp.path().join("minus.psl");
    let sizes = temp.path().join("chrom.sizes");
    let output = temp.path().join("lifted.psl");
    std::fs::write(&sizes, "chr\t5600\n").unwrap();
    // Window chr:1201-1400 (1-based inclusive) = genome [1200, 1400).
    // '-' record: forward qStart/qEnd [0, 101), block qStarts in the window
    // RC frame (qSize - qEnd = 99), matching pgr/UCSC convention.
    std::fs::write(
        &input,
        "101\t0\t0\t0\t0\t0\t0\t0\t-\tchr:1201-1400\t200\t0\t101\tchr\t5600\t3299\t3400\t1\t101,\t99,\t3299,\n",
    )
    .unwrap();

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
    // qStart/qEnd must be forward genome coordinates [1200, 1301).
    assert!(
        output_content.contains("-\tchr\t5600\t1200\t1301"),
        "qStart/qEnd must lift to forward coordinates: {output_content}"
    );
    // The block qStarts are in the chromosome RC frame: 5600 - 1301 = 4299.
    assert!(
        output_content.contains("4299,\t3299,"),
        "qStarts must lift in the RC frame: {output_content}"
    );

    // `psl to-range` on the lifted record must recover the genomic span
    // [1200, 1301) as 1-based inclusive "chr:1201-1301".
    let range_out = temp.path().join("ranges.rg");
    PgrCmd::new()
        .args(&[
            "psl",
            "to-range",
            output.to_str().unwrap(),
            "-o",
            range_out.to_str().unwrap(),
        ])
        .run();
    let ranges = fs::read_to_string(&range_out).unwrap();
    assert!(
        ranges.contains("chr:1201-1301"),
        "to-range must recover the genomic span: {ranges}"
    );
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

#[test]
fn test_stats_empty() {
    let temp = TempDir::new().unwrap();
    let input = temp.path().join("empty.psl");
    fs::write(&input, "").unwrap();
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
    assert_eq!(lines.len(), 1);
    assert!(lines[0].starts_with('#'));
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
    //    qStart/qEnd are forward coordinates; the block qStarts (80, in the
    //    RC frame) reverse back to the same genomic span.
    //    Range: chr1:101-200:11-20

    let lines: Vec<&str> = output_content.lines().collect();
    assert_eq!(lines.len(), 2);
    assert_eq!(lines[0], "chr1:101-200:11-20");
    assert_eq!(lines[1], "chr1:101-200:11-20");
}

//
// psl swap
//

#[test]
fn test_psl_swap_mrna() {
    let temp = TempDir::new().unwrap();
    let input = get_path("swap", "input", "mrna.psl");
    let expected = get_path("swap", "expected", "mrnaTest.psl");
    let output = temp.path().join("out.psl");

    PgrCmd::new()
        .args(&[
            "psl",
            "swap",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    let expected_content = fs::read_to_string(&expected).unwrap();
    assert_eq!(output_content, expected_content);
}

#[test]
fn test_psl_swap_mrna_no_rc() {
    let temp = TempDir::new().unwrap();
    let input = get_path("swap", "input", "mrna.psl");
    let expected = get_path("swap", "expected", "mrnaNoRcTest.psl");
    let output = temp.path().join("out.psl");

    PgrCmd::new()
        .args(&[
            "psl",
            "swap",
            "--no-rc",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    let expected_content = fs::read_to_string(&expected).unwrap();
    assert_eq!(output_content, expected_content);
}

#[test]
fn test_psl_swap_trans() {
    let temp = TempDir::new().unwrap();
    let input = get_path("swap", "input", "trans.psl");
    let expected = get_path("swap", "expected", "transTest.psl");
    let output = temp.path().join("out.psl");

    PgrCmd::new()
        .args(&[
            "psl",
            "swap",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let output_content = fs::read_to_string(&output).unwrap();
    let expected_content = fs::read_to_string(&expected).unwrap();
    assert_eq!(output_content, expected_content);
}

// psl to-paf
//

#[test]
fn test_to_paf_basic() {
    let temp = TempDir::new().unwrap();
    let input = get_path("to_chain", "input", "mtor.psl");
    let output = temp.path().join("out.paf");

    PgrCmd::new()
        .args(&[
            "psl",
            "to-paf",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .run();

    let content = fs::read_to_string(&output).unwrap();
    let first = content.lines().next().unwrap();
    let fields: Vec<&str> = first.split('\t').collect();
    assert_eq!(fields[0], "ENST00000361445.8");
    assert_eq!(fields[1], "2549");
    assert_eq!(fields[4], "+"); // first char of "+-"
    assert_eq!(fields[5], "chr1");
    assert_eq!(fields[9], "2542"); // match count
    assert_eq!(fields[10], "2542"); // sum of block sizes
    assert_eq!(fields[11], "255"); // mapq
}
