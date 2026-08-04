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
fn command_rg_count() {
    let dir = TempDir::new().unwrap();
    let target = dir.path().join("target.rg");
    let intervals = dir.path().join("iv.rg");
    std::fs::write(&target, "chr1:1-10\nchr1:5-15\nchr2:100-200\nbad line\n").unwrap();
    std::fs::write(&intervals, "chr1:1-10\nchr1:5-15\nchr2:150-160\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "rg",
            "count",
            target.to_str().unwrap(),
            intervals.to_str().unwrap(),
        ])
        .run();
    assert_eq!(stdout, "chr1:1-10\t2\nchr1:5-15\t2\nchr2:100-200\t1\n");
}

#[test]
fn command_rg_count_stdin_intervals() {
    let dir = TempDir::new().unwrap();
    let target = dir.path().join("target.rg");
    std::fs::write(&target, "chr1:1-5\nchr1:10-20\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "count", target.to_str().unwrap(), "stdin"])
        .stdin("chr1:1-10\nchr1:15-16\n")
        .run();
    assert_eq!(stdout, "chr1:1-5\t1\nchr1:10-20\t2\n");
}

// Migrated from intspan `tests/cli_rgr.rs` `command_count` (.rg part; the
// TSV `-f` part tests functionality dropped by the rg-family input contract).
#[test]
fn command_rg_count_fixture() {
    let (stdout, _) = cmd(&["count", &fixture("S288c.rg"), &fixture("S288c.rg")]).run();
    assert_eq!(stdout.lines().count(), 6);
    assert!(stdout.contains("I:1-100\t2"), "got: {stdout}");
    assert!(stdout.contains("21294-22075\t1"), "got: {stdout}");
}

#[test]
fn command_rg_prop() {
    let dir = TempDir::new().unwrap();
    let json = dir.path().join("in.json");
    let rg = dir.path().join("a.rg");
    std::fs::write(&json, r#"{"chr1":"1-10,20-30"}"#).unwrap();
    std::fs::write(&rg, "chr1:5-25\nchr1:50-60\nchr2:1-10\nbad line\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "prop", json.to_str().unwrap(), rg.to_str().unwrap()])
        .run();
    // chr1:5-25 intersects 5-10 (6 bp) + 20-25 (6 bp) = 12 of 21 bp.
    assert_eq!(
        stdout,
        "chr1:5-25\t0.5714\nchr1:50-60\t0.0000\nchr2:1-10\t0.0000\n"
    );
}

#[test]
fn command_rg_prop_full() {
    let dir = TempDir::new().unwrap();
    let json = dir.path().join("in.json");
    let rg = dir.path().join("a.rg");
    std::fs::write(&json, r#"{"chr1":"1-10"}"#).unwrap();
    std::fs::write(&rg, "chr1:1-20\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "rg",
            "prop",
            json.to_str().unwrap(),
            rg.to_str().unwrap(),
            "--full",
        ])
        .run();
    assert_eq!(stdout, "chr1:1-20\t0.5000\t20\t10\n");
}

// Migrated from intspan `tests/cli_rgr.rs` `command_prop` (.rg part).
#[test]
fn command_rg_prop_fixture() {
    let (stdout, _) = cmd(&["prop", &fixture("intergenic.json"), &fixture("S288c.rg")]).run();
    assert_eq!(stdout.lines().count(), 6);
    assert!(stdout.contains("I:1-100\t0.0000"), "got: {stdout}");
    assert!(stdout.contains("II:21294-22075\t1.0000"), "got: {stdout}");
}

#[test]
fn command_rg_sort() {
    let dir = TempDir::new().unwrap();
    let rg = dir.path().join("a.rg");
    std::fs::write(
        &rg,
        "chr2:100-200\nchr1:50-60\nchr1:100-110\nbad line\nchr1:20-30\n",
    )
    .unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "sort", rg.to_str().unwrap()])
        .run();
    assert_eq!(
        stdout,
        "chr1:20-30\nchr1:50-60\nchr1:100-110\nchr2:100-200\nbad line\n"
    );
}

#[test]
fn command_rg_sort_strand_order() {
    let dir = TempDir::new().unwrap();
    let rg = dir.path().join("a.rg");
    std::fs::write(&rg, "chr1(+):50-60\nchr1(-):50-60\nchr1(+):20-30\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["rg", "sort", rg.to_str().unwrap()])
        .run();
    // Key is (chr, start, strand): start 20 first, then start 50 with "+"
    // before "-" (ASCII '+' < '-').
    assert_eq!(stdout, "chr1(+):20-30\nchr1(+):50-60\nchr1(-):50-60\n");
}

// Migrated from intspan `tests/cli_rgr.rs` `command_sort` (.rg part).
#[test]
fn command_rg_sort_fixture() {
    let (stdout, _) = cmd(&["sort", &fixture("S288c.rg")]).run();
    assert_eq!(stdout.lines().count(), 6);
    assert_eq!(stdout.lines().next().unwrap(), "I:1-100");
    assert_eq!(stdout.lines().last().unwrap(), "II:23537-24097");
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
