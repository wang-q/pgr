//! Migrated from the external intspan project's `tests/cli_spanr.rs`
//! (spanr CLI), adapted to run against `pgr runlist`. Fixtures live in
//! `tests/runlist/` (copied from `intspan/tests/spanr/`).

#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::path::PathBuf;

fn fixture(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/runlist")
        .join(name)
        .to_string_lossy()
        .into_owned()
}

fn cmd(args: &[&str]) -> PgrCmd {
    let mut full = vec!["runlist"];
    full.extend_from_slice(args);
    PgrCmd::new().args(&full)
}

#[test]
fn command_invalid() {
    let (_, stderr) = PgrCmd::new().args(&["runlist", "foobar"]).run_fail();
    assert!(stderr.contains("recognized"), "got: {stderr}");
}

#[test]
fn bare_runlist_shows_help() {
    let (_, stderr) = PgrCmd::new().args(&["runlist"]).run_fail();
    assert!(stderr.contains("Usage: pgr runlist"), "got: {stderr}");
}

#[test]
fn file_doesnt_provided() {
    let (_, stderr) = cmd(&["genome"]).run_fail();
    assert!(stderr.contains("not provided"), "got: {stderr}");
}

#[test]
fn file_doesnt_exist() {
    let (_, stderr) = cmd(&["genome", "tests/file/doesnt/exist"]).run_fail();
    assert!(stderr.contains("could not open"), "got: {stderr}");
}

#[test]
fn command_genome() {
    let (stdout, _) = cmd(&["genome", &fixture("S288c.chr.sizes")]).run();
    assert!(stdout.contains("1-230218"), "got: {stdout}");
}

#[test]
fn command_some() {
    let (stdout, _) = cmd(&["some", &fixture("Atha.json"), &fixture("Atha.list")]).run();
    let lines = stdout.lines().count();
    assert!(lines == 11 || lines == 12, "line count {lines}");
    assert!(stdout.contains("AT2G01008"), "got: {stdout}");
    assert!(!stdout.contains("AT2G01021"), "got: {stdout}");
}

#[test]
fn command_merge() {
    let (stdout, _) = cmd(&["merge", &fixture("I.json"), &fixture("II.json")]).run();
    let lines = stdout.lines().count();
    assert!(lines == 8 || lines == 9, "line count {lines}");
    assert!(stdout.contains("28547-29194"), "got: {stdout}");
    assert!(stdout.contains("\"II\":"), "got: {stdout}");

    let (stdout, _) = cmd(&[
        "merge",
        &fixture("I.json"),
        &fixture("II.other.json"),
        "--all",
    ])
    .run();
    let lines = stdout.lines().count();
    assert!(lines == 8 || lines == 9, "line count {lines}");
    assert!(stdout.contains("28547-29194"), "got: {stdout}");
    assert!(stdout.contains("\"II.other\":"), "got: {stdout}");
}

#[test]
fn command_split() {
    let (stdout, _) = cmd(&["split", &fixture("I.II.json")]).run();
    assert!(stdout.contains("28547-29194"), "got: {stdout}");
    assert!(stdout.contains("{\"I\":"), "got: {stdout}");
    assert!(stdout.contains("{\"II\":"), "got: {stdout}");
}

#[test]
fn command_split_to() {
    let dir = tempfile::TempDir::new().unwrap();
    let (stdout, _) = cmd(&[
        "split",
        &fixture("I.II.json"),
        "-o",
        dir.path().to_str().unwrap(),
    ])
    .run();
    assert_eq!(stdout, "");
    assert!(dir.path().join("II.json").is_file());
    assert!(!dir.path().join("I.II.json").exists());
}

#[test]
fn command_stat() {
    let (stdout, _) = cmd(&[
        "stat",
        &fixture("S288c.chr.sizes"),
        &fixture("intergenic.json"),
    ])
    .run();
    assert_eq!(stdout.lines().count(), 18, "line count");
    assert_eq!(
        stdout.lines().next().unwrap().split(',').count(),
        4,
        "field count"
    );
    assert!(stdout.contains("all,12071326,1059702,"), "got: {stdout}");
}

#[test]
fn command_stat_all() {
    let (stdout, _) = cmd(&[
        "stat",
        &fixture("S288c.chr.sizes"),
        &fixture("intergenic.json"),
        "--all",
    ])
    .run();
    assert_eq!(stdout.lines().count(), 2, "line count");
    assert_eq!(
        stdout.lines().next().unwrap().split(',').count(),
        3,
        "field count"
    );
    assert!(!stdout.contains("all"), "got: {stdout}");
}

#[test]
fn command_statop() {
    let (stdout, _) = cmd(&[
        "statop",
        &fixture("S288c.chr.sizes"),
        &fixture("intergenic.json"),
        &fixture("repeat.json"),
    ])
    .run();
    assert_eq!(stdout.lines().count(), 18, "line count");
    assert_eq!(
        stdout.lines().next().unwrap().split(',').count(),
        8,
        "field count"
    );
    assert!(stdout.contains("36721"), "sum exists: {stdout}");
    assert!(stdout.contains(",repeatLength,"), "got: {stdout}");
    assert!(stdout.contains("\nI,"), "got: {stdout}");
    assert!(stdout.contains("\nXVI,"), "got: {stdout}");
}

#[test]
fn command_statop_all() {
    let (stdout, _) = cmd(&[
        "statop",
        &fixture("S288c.chr.sizes"),
        &fixture("intergenic.json"),
        &fixture("repeat.json"),
        "--all",
    ])
    .run();
    assert_eq!(stdout.lines().count(), 2, "line count");
    assert_eq!(
        stdout.lines().next().unwrap().split(',').count(),
        7,
        "field count"
    );
    assert!(stdout.contains("36721"), "sum exists: {stdout}");
    assert!(stdout.contains(",repeatLength,"), "got: {stdout}");
    assert!(!stdout.contains("\nI,"), "got: {stdout}");
    assert!(!stdout.contains("\nXVI,"), "got: {stdout}");
}

#[test]
fn command_statop_invalid() {
    let (_, stderr) = cmd(&[
        "statop",
        &fixture("S288c.chr.sizes"),
        &fixture("intergenic.json"),
        &fixture("repeat.json"),
        "--op",
        "invalid",
        "--all",
    ])
    .run_fail();
    assert!(stderr.contains("invalid value"), "got: {stderr}");
}

#[test]
fn command_combine() {
    let (stdout, _) = cmd(&["combine", &fixture("Atha.json")]).run();
    let lines = stdout.lines().count();
    assert!(lines == 4 || lines == 5, "line count {lines}");
    assert!(!stdout.contains("7232,7384"), "combined: {stdout}");

    let (stdout, _) = cmd(&["combine", &fixture("Atha.json"), "--op", "xor"]).run();
    let lines = stdout.lines().count();
    assert!(lines == 4 || lines == 5, "line count {lines}");
    assert!(stdout.contains("7233-7383"), "xor: {stdout}");

    let (stdout, _) = cmd(&["combine", &fixture("II.json")]).run();
    let lines = stdout.lines().count();
    assert!(lines == 2 || lines == 3, "line count {lines}");
    assert!(stdout.contains("21294-22075,"), "no changes: {stdout}");
}

#[test]
fn command_compare() {
    let (stdout, _) = cmd(&[
        "compare",
        &fixture("intergenic.json"),
        &fixture("repeat.json"),
        "--op",
        "intersect",
    ])
    .run();
    let lines = stdout.lines().count();
    assert!(lines == 18 || lines == 19, "line count {lines}");
    assert!(stdout.contains("878539-878709"), "runlist exists: {stdout}");
    assert!(stdout.contains("\"XVI\":"), "got: {stdout}");

    let (stdout, _) = cmd(&[
        "compare",
        &fixture("intergenic.json"),
        &fixture("repeat.json"),
        "--op",
        "union",
    ])
    .run();
    let lines = stdout.lines().count();
    assert!(lines == 18 || lines == 19, "line count {lines}");
    assert!(!stdout.contains("\"-\""), "no empty runlists: {stdout}");
    assert!(stdout.contains("\"XVI\":"), "got: {stdout}");

    let (stdout, _) = cmd(&[
        "compare",
        &fixture("intergenic.json"),
        &fixture("repeat.json"),
        "--op",
        "xor",
    ])
    .run();
    let lines = stdout.lines().count();
    assert!(lines == 18 || lines == 19, "line count {lines}");
    assert!(!stdout.contains("\"-\""), "no empty runlists: {stdout}");
    assert!(stdout.contains("\"XVI\":"), "got: {stdout}");

    let (stdout, _) = cmd(&[
        "compare",
        &fixture("I.II.json"),
        &fixture("repeat.json"),
        "--op",
        "intersect",
    ])
    .run();
    let lines = stdout.lines().count();
    assert!(lines == 38 || lines == 39, "line count {lines}");

    let (stdout, _) = cmd(&[
        "compare",
        &fixture("I.II.json"),
        &fixture("I.json"),
        &fixture("II.json"),
        "--op",
        "intersect",
    ])
    .run();
    let lines = stdout.lines().count();
    assert!(lines == 10 || lines == 11, "line count {lines}");
    assert!(!stdout.contains("13744-17133"), "all empty: {stdout}");
}

#[test]
fn command_span() {
    let (stdout, _) = cmd(&["span", &fixture("brca2.json"), "--op", "cover"]).run();
    let lines = stdout.lines().count();
    assert!(lines == 2 || lines == 3, "line count {lines}");
    assert!(stdout.contains("32316461-32398770"), "cover: {stdout}");

    let (stdout, _) = cmd(&["span", &fixture("brca2.json"), "--op", "fill", "-n", "1000"]).run();
    let lines = stdout.lines().count();
    assert!(lines == 2 || lines == 3, "line count {lines}");
    assert!(
        stdout.contains("32325076-32326613"),
        "newly emerged: {stdout}"
    );
    assert_ne!(stdout.matches(',').count(), 25, "original");
    assert_eq!(stdout.matches(',').count(), 18, "new");

    let (stdout, _) = cmd(&["span", &fixture("brca2.json"), "--op", "trim", "-n", "200"]).run();
    let lines = stdout.lines().count();
    assert!(lines == 2 || lines == 3, "line count {lines}");
    assert_ne!(stdout.matches(',').count(), 25, "original");
    assert_eq!(stdout.matches(',').count(), 3, "new");

    let (stdout, _) = cmd(&["span", &fixture("brca2.json"), "--op", "pad", "-n", "2000"]).run();
    let lines = stdout.lines().count();
    assert!(lines == 2 || lines == 3, "line count {lines}");
    assert_ne!(stdout.matches(',').count(), 25, "original");
    assert_eq!(stdout.matches(',').count(), 6, "new");

    let (stdout, _) = cmd(&[
        "span",
        &fixture("brca2.json"),
        "--op",
        "excise",
        "-n",
        "400",
    ])
    .run();
    let lines = stdout.lines().count();
    assert!(lines == 2 || lines == 3, "line count {lines}");
    assert_ne!(stdout.matches(',').count(), 25, "original");
    assert_eq!(stdout.matches(',').count(), 3, "new");
}

#[test]
fn command_span_invalid() {
    let (_, stderr) = cmd(&["span", &fixture("brca2.json"), "--op", "invalid"]).run_fail();
    assert!(stderr.contains("invalid value"), "got: {stderr}");
}

#[test]
fn command_convert() {
    let (stdout, _) = cmd(&["convert", &fixture("repeat.json")]).run();
    assert_eq!(stdout.lines().count(), 28, "line count");
    assert!(stdout.contains("II:327069-327703"), "first range: {stdout}");

    let (stdout, _) = cmd(&["convert", &fixture("repeat.json"), "--longest"]).run();
    assert_eq!(stdout.lines().count(), 11, "line count");
    assert!(stdout.contains("IV:981142-987119"), "longest: {stdout}");
    assert!(
        !stdout.contains("IV:757572-759779"),
        "not longest: {stdout}"
    );
}
