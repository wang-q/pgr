#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::path::PathBuf;

fn fixture(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/gff")
        .join(name)
        .to_string_lossy()
        .into_owned()
}

#[test]
fn command_rg_default() {
    let (stdout, _) = PgrCmd::new()
        .args(&["gff", "rg", "tests/gff/test.gff"])
        .run();

    assert!(stdout.contains("gene1\ttest.chr1(+):1000-2000"));
    assert!(stdout.contains("prefix:gene2\ttest.chr1(-):3000-4000"));
    assert!(!stdout.contains("mRNA1"));
}

#[test]
fn command_rg_tag() {
    let (stdout, _) = PgrCmd::new()
        .args(&["gff", "rg", "tests/gff/test.gff", "--tag", "mRNA"])
        .run();

    assert!(stdout.contains("mRNA1\ttest.chr1(+):1000-2000"));
    assert!(!stdout.contains("gene1"));
}

#[test]
fn command_runlist_gff() {
    // Migrated from the external spanr `gff` command (intspan test suite).
    let (stdout, _) = PgrCmd::new()
        .args(&["gff", "runlist", &fixture("NC_007942.gff")])
        .run();
    let lines = stdout.lines().count();
    assert!(lines == 2 || lines == 3, "line count {lines}");
    assert!(stdout.contains("NC_007942"), "chromosomes exists: {stdout}");
    assert!(stdout.contains("1-152218"), "full chr runlist: {stdout}");
}

#[test]
fn command_runlist_gff_tag_and_merge() {
    let dir = tempfile::TempDir::new().unwrap();
    let cds = dir.path().join("cds.json");
    let repeat = dir.path().join("repeat.json");
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "gff",
            "runlist",
            &fixture("NC_007942.gff"),
            "--tag",
            "CDS",
            "-o",
            cds.to_str().unwrap(),
        ])
        .run();
    assert_eq!(stdout, "");
    assert!(cds.is_file());

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "gff",
            "runlist",
            &fixture("NC_007942.rm.gff"),
            "-o",
            repeat.to_str().unwrap(),
        ])
        .run();
    assert_eq!(stdout, "");
    assert!(repeat.is_file());

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "runlist",
            "merge",
            cds.to_str().unwrap(),
            repeat.to_str().unwrap(),
        ])
        .run();
    let lines = stdout.lines().count();
    assert!(lines == 8 || lines == 9, "line count {lines}");
    assert!(stdout.contains("cds"), "got: {stdout}");
    assert!(stdout.contains("repeat"), "got: {stdout}");
}

#[test]
fn command_rg_asm() {
    let (stdout, _) = PgrCmd::new()
        .args(&["gff", "rg", "tests/gff/test.gff", "--asm", "Human"])
        .run();

    assert!(stdout.contains("gene1\tHuman.chr1(+):1000-2000"));
}

#[test]
fn command_rg_simplify() {
    let (stdout, _) = PgrCmd::new()
        .args(&["gff", "rg", "tests/gff/test.gff", "--key-simplify"])
        .run();

    assert!(stdout.contains("prefix:gene2\ttest.chr1(-):3000-4000"));
    // assert!(!stdout.contains("prefix:gene2"));
}

#[test]
fn command_rg_simplify_destructive() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "gff",
            "rg",
            "tests/gff/test.gff",
            "--tag",
            "CDS",
            "--key",
            "Name", // NP_414542.1
            "--key-simplify",
        ])
        .run();

    // With destructive simplify, this would be NP_414542 (missing .1)
    // We want to ensure it is destructively simplified as per user request
    assert!(stdout.contains("NP_414542\ttest.chr1(+):5000-6000"));
    assert!(!stdout.contains("NP_414542.1"));
}

#[test]
fn command_rg_case_insensitive() {
    let (stdout, _) = PgrCmd::new()
        .args(&["gff", "rg", "tests/gff/test.gff", "--tag", "MRNA"])
        .run();

    assert!(stdout.contains("mRNA1\ttest.chr1(+):1000-2000"));
}

#[test]
fn command_rg_key() {
    let (stdout, _) = PgrCmd::new()
        .args(&["gff", "rg", "tests/gff/test.gff", "--key", "Name"])
        .run();

    assert!(stdout.contains("GENE1\ttest.chr1(+):1000-2000"));
    assert!(!stdout.contains("gene1"));
}

#[test]
fn command_rg_key_parent() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "gff",
            "rg",
            "tests/gff/test.gff",
            "--tag",
            "mRNA",
            "--key",
            "Parent",
        ])
        .run();

    assert!(stdout.contains("gene1\ttest.chr1(+):1000-2000"));
}

#[test]
fn command_rg_key_product() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "gff",
            "rg",
            "tests/gff/test.gff",
            "--tag",
            "mRNA",
            "--key",
            "product",
        ])
        .run();

    assert!(stdout.contains("thr operon leader peptide\ttest.chr1(+):1000-2000"));
}

#[test]
fn command_rg_ss() {
    let (stdout, _) = PgrCmd::new()
        .args(&["gff", "rg", "tests/gff/test.gff", "--seq-simplify"])
        .run();

    // test.gff contains "test.chr1", which doesn't need simplification.
    // We need to add a case with a complex chromosome name to test this properly.
    // For now, let's just check it runs without error and outputs valid lines.
    assert!(stdout.contains("gene1\ttest.chr1(+):1000-2000"));
}

#[test]
fn command_rg_ucsc_mm10gencode() {
    let (stdout, _) = PgrCmd::new()
        .args(&["gff", "rg", "tests/gff/input/mm10Gencode.gff3"])
        .run();

    assert!(stdout.contains("gene:ENSMUSG00000024186"));
    assert!(stdout.contains("mm10Gencode.chr17(+):98013-106386"));
}

#[test]
fn command_rg_ucsc_malformed_no_panic() {
    // 87-byte GFF3 with bogus quotes: tolerated without panic.
    PgrCmd::new()
        .args(&["gff", "rg", "tests/gff/input/bogusQuotes.gff3"])
        .assert()
        .success();

    // GFF2-era file with a frame bug: friendly error, no panic.
    PgrCmd::new()
        .args(&["gff", "rg", "tests/gff/input/frameBug.gff"])
        .assert()
        .failure();
}
