#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

// ── paf validate (CIGAR self-consistency) ─────────────

const VALID: &str = "q\t100\t0\t50\t+\tt\t200\t0\t50\t45\t50\t255\tcg:Z:50M\n";

#[test]
fn validate_clean_paf_reports_no_invalid() {
    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "validate", "stdin"])
        .stdin(VALID)
        .run();
    assert!(stderr.is_empty(), "stderr: {stderr}");
    assert!(stdout.contains("Total records: 1"), "stdout: {stdout}");
    assert!(stdout.contains("Query invalid records: 0"));
    assert!(stdout.contains("Target invalid records: 0"));
}

#[test]
fn validate_flags_bad_query_and_target_ends() {
    // query_end (40) < CIGAR-derived query span, target_end (60) > derived.
    let paf = "q\t100\t0\t40\t+\tt\t200\t0\t60\t10\t50\t255\tcg:Z:10M5I5D\n";
    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "validate", "stdin"])
        .stdin(paf)
        .run();
    assert!(stderr.is_empty(), "stderr: {stderr}");
    assert!(stdout.contains("Query invalid records: 1"));
    assert!(stdout.contains("Target invalid records: 1"));
    assert!(stdout.contains("q:0-40"));
    assert!(stdout.contains("t:0-60"));
}

#[test]
fn validate_counts_records_without_cigar() {
    // No cg:Z tag -> counted, not fatal.
    let paf = "q\t100\t0\t50\t+\tt\t200\t0\t50\t45\t50\t255\n";
    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "validate", "stdin"])
        .stdin(paf)
        .run();
    assert!(stderr.is_empty(), "stderr: {stderr}");
    assert!(stdout.contains("Records without cg:Z tag: 1"));
    assert!(stdout.contains("Query invalid records: 0"));
}

#[test]
fn validate_handles_malformed_cigar_without_panicking() {
    // 'N' is not a valid PAF CIGAR op -> counted as malformed, not fatal.
    let paf = "q\t100\t0\t50\t+\tt\t200\t0\t50\t45\t50\t255\tcg:Z:50N\n";
    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "validate", "stdin"])
        .stdin(paf)
        .run();
    assert!(stderr.is_empty(), "stderr: {stderr}");
    assert!(stdout.contains("Records with malformed cg:Z tag: 1"));
}

#[test]
fn validate_skips_comments_and_empty_lines() {
    let paf = format!("# header\n\n{VALID}\n# comment\n");
    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "validate", "stdin"])
        .stdin(paf)
        .run();
    assert!(stderr.is_empty(), "stderr: {stderr}");
    assert!(stdout.contains("Total records: 1"));
}

#[test]
fn validate_writes_report_to_outfile() {
    let dir = std::env::temp_dir().join(format!("pgr_validate_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let out = dir.join("report.txt");
    let (stdout, stderr) = PgrCmd::new()
        .args(&["paf", "validate", "stdin", "-o", out.to_str().unwrap()])
        .stdin(VALID)
        .run();
    assert!(stderr.is_empty(), "stderr: {stderr}");
    assert!(stdout.is_empty(), "stdout should be empty: {stdout}");
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.contains("Total records: 1"));
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn validate_rejects_outfile_overwriting_input() {
    let dir = std::env::temp_dir().join(format!("pgr_validate_in_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let paf_path = dir.join("in.paf");
    std::fs::write(&paf_path, VALID).unwrap();
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "validate",
            paf_path.to_str().unwrap(),
            "-o",
            paf_path.to_str().unwrap(),
        ])
        .run_fail();
    assert!(stderr.contains("is also an input file"), "stderr: {stderr}");
    std::fs::remove_dir_all(&dir).ok();
}
