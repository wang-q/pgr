//! CLI tests for `pgr rg` (migrated from `pgr runlist cover/coverage`).

#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::path::PathBuf;
use tempfile::TempDir;

fn fixture_dir() -> (TempDir, PathBuf) {
    let dir = TempDir::new().unwrap();
    let rg = dir.path().join("a.rg");
    std::fs::write(&rg, "chr1:1-10\nchr1:5-15\nchr2(+):100-200\nbad line\n").unwrap();
    let rg2 = dir.path().join("b.rg");
    std::fs::write(&rg2, "chr1:20-25\nchr2:150-160\n").unwrap();
    (dir, rg)
}

fn fixture(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/runlist")
        .join(name)
        .to_string_lossy()
        .into_owned()
}

fn cmd(args: &[&str]) -> PgrCmd {
    let mut full = vec!["rg"];
    full.extend_from_slice(args);
    PgrCmd::new().args(&full)
}

#[test]
fn command_rg_cover() {
    let (dir, rg) = fixture_dir();
    let rg2 = dir.path().join("b.rg");
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "cover", rg.to_str().unwrap(), rg2.to_str().unwrap()])
        .run();
    assert_eq!(
        stdout,
        "{\n  \"chr1\": \"1-15,20-25\",\n  \"chr2\": \"100-200\"\n}\n"
    );
}

#[test]
fn command_rg_cover_stdin() {
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "cover", "stdin"])
        .stdin("chr1:1-5\nchr1:3-8\n")
        .run();
    assert_eq!(stdout, "{\n  \"chr1\": \"1-8\"\n}\n");
}

#[test]
fn command_rg_cover_skips_reversed_ranges() {
    // `chr1:10-5` (start > end) must be skipped, not panic.
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "cover", "stdin"])
        .stdin("chr1:10-5\nchr1:1-10\n")
        .run();
    assert_eq!(stdout, "{\n  \"chr1\": \"1-10\"\n}\n");
    // Coordinates above the representable maximum must be skipped too.
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "cover", "stdin"])
        .stdin("chr1:2147483647-2147483647\n")
        .run();
    assert_eq!(stdout, "{}\n");
}

#[test]
fn command_rg_coverage() {
    let (dir, rg) = fixture_dir();
    let rg2 = dir.path().join("b.rg");
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "rg",
            "coverage",
            rg.to_str().unwrap(),
            rg2.to_str().unwrap(),
            "-m",
            "2",
        ])
        .run();
    assert_eq!(
        stdout,
        "{\n  \"chr1\": \"5-10\",\n  \"chr2\": \"150-160\"\n}\n"
    );
}

#[test]
fn command_rg_coverage_detailed() {
    let (_dir, rg) = fixture_dir();
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "coverage", rg.to_str().unwrap(), "-m", "1", "-d"])
        .run();
    assert_eq!(
        stdout,
        "{\n  \"1\": {\n    \"chr1\": \"1-4,11-15\",\n    \"chr2\": \"100-200\"\n  },\n  \"2\": {\n    \"chr1\": \"5-10\"\n  }\n}\n"
    );
}

#[test]
fn command_cover() {
    let (stdout, _) = cmd(&["cover", &fixture("S288c.rg")]).run();
    let lines = stdout.lines().count();
    assert!(lines == 3 || lines == 4, "line count {lines}");
    assert!(!stdout.contains("S288c"), "species name: {stdout}");
    assert!(!stdout.contains("1-100"), "merged: {stdout}");
    assert!(stdout.contains("1-150"), "covered: {stdout}");

    let (stdout, _) = cmd(&["cover", &fixture("dazzname.rg")]).run();
    let lines = stdout.lines().count();
    assert!(lines == 2 || lines == 3, "line count {lines}");
    assert!(stdout.contains("infile_0/1/0_514"), "chr name: {stdout}");
    assert!(stdout.contains("19-499"), "covered: {stdout}");
}

#[test]
fn command_coverage() {
    let (stdout, _) = cmd(&["coverage", &fixture("S288c.rg"), "-m", "2"]).run();
    let lines = stdout.lines().count();
    assert!(lines == 3 || lines == 4, "line count {lines}");
    assert!(!stdout.contains("S288c"), "species name: {stdout}");
    assert!(!stdout.contains("1-150"), "coverage 1: {stdout}");
    assert!(stdout.contains("90-100"), "coverage 2: {stdout}");
}

#[test]
fn command_coverage_detailed() {
    let (stdout, _) = cmd(&["coverage", &fixture("S288c.rg"), "-m", "1", "-d"]).run();
    let lines = stdout.lines().count();
    assert!(lines == 9 || lines == 10, "line count {lines}");
    assert!(!stdout.contains("S288c"), "species name: {stdout}");
    assert!(stdout.contains("1-89"), "coverage 1: {stdout}");
    assert!(stdout.contains("90-100"), "coverage 2: {stdout}");
    assert!(stdout.contains("190-200"), "coverage 2: {stdout}");
}
