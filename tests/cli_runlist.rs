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

#[test]
fn command_runlist_cover() {
    let (dir, rg) = fixture_dir();
    let rg2 = dir.path().join("b.rg");
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "runlist",
            "cover",
            rg.to_str().unwrap(),
            rg2.to_str().unwrap(),
        ])
        .run();
    assert_eq!(
        stdout,
        "{\n  \"chr1\": \"1-15,20-25\",\n  \"chr2\": \"100-200\"\n}\n"
    );
}

#[test]
fn command_runlist_cover_stdin() {
    let (stdout, _) = PgrCmd::new()
        .args(&["runlist", "cover", "stdin"])
        .stdin("chr1:1-5\nchr1:3-8\n")
        .run();
    assert_eq!(stdout, "{\n  \"chr1\": \"1-8\"\n}\n");
}

#[test]
fn command_runlist_coverage() {
    let (dir, rg) = fixture_dir();
    let rg2 = dir.path().join("b.rg");
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "runlist",
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
fn command_runlist_coverage_detailed() {
    let (_dir, rg) = fixture_dir();
    let (stdout, _) = PgrCmd::new()
        .args(&["runlist", "coverage", rg.to_str().unwrap(), "-m", "1", "-d"])
        .run();
    assert_eq!(
        stdout,
        "{\n  \"1\": {\n    \"chr1\": \"1-4,11-15\",\n    \"chr2\": \"100-200\"\n  },\n  \"2\": {\n    \"chr1\": \"5-10\"\n  }\n}\n"
    );
}

#[test]
fn command_runlist_span_fill_excise() {
    let (dir, _rg) = fixture_dir();
    let json = dir.path().join("in.json");
    std::fs::write(&json, r#"{"chr1":"1-3,7-10,15-16"}"#).unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "runlist",
            "span",
            json.to_str().unwrap(),
            "--op",
            "fill",
            "-n",
            "3",
        ])
        .run();
    assert_eq!(stdout, "{\n  \"chr1\": \"1-10,15-16\"\n}\n");
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "runlist",
            "span",
            json.to_str().unwrap(),
            "--op",
            "excise",
            "-n",
            "3",
        ])
        .run();
    assert_eq!(stdout, "{\n  \"chr1\": \"1-3,7-10\"\n}\n");
}

#[test]
fn command_runlist_compare_intersect() {
    let dir = TempDir::new().unwrap();
    let a = dir.path().join("a.json");
    let b = dir.path().join("b.json");
    std::fs::write(&a, r#"{"chr1":"1-15,20-25","chr2":"100-200"}"#).unwrap();
    std::fs::write(&b, r#"{"chr1":"5-10","chr2":"150-160"}"#).unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "runlist",
            "compare",
            a.to_str().unwrap(),
            b.to_str().unwrap(),
            "--op",
            "intersect",
        ])
        .run();
    assert_eq!(
        stdout,
        "{\n  \"chr1\": \"5-10\",\n  \"chr2\": \"150-160\"\n}\n"
    );
}

#[test]
fn command_runlist_merge() {
    let dir = TempDir::new().unwrap();
    let a = dir.path().join("sample.a.json");
    let b = dir.path().join("sample.b.json");
    std::fs::write(&a, r#"{"chr1":"1-5"}"#).unwrap();
    std::fs::write(&b, r#"{"chr2":"6-9"}"#).unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["runlist", "merge", a.to_str().unwrap(), b.to_str().unwrap()])
        .run();
    // Both stems start with "sample": only one key survives (spanr parity).
    assert_eq!(
        stdout,
        "{\n  \"sample\": {\n    \"chr2\": \"6-9\"\n  }\n}\n"
    );
}

#[test]
fn command_runlist_genome() {
    let dir = TempDir::new().unwrap();
    let sizes = dir.path().join("sizes.txt");
    std::fs::write(&sizes, "chr1\t1000\nchr2\t2000\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["runlist", "genome", sizes.to_str().unwrap()])
        .run();
    assert_eq!(
        stdout,
        "{\n  \"chr1\": \"1-1000\",\n  \"chr2\": \"1-2000\"\n}\n"
    );
}

#[test]
fn command_runlist_combine() {
    let dir = TempDir::new().unwrap();
    let json = dir.path().join("multi.json");
    std::fs::write(
        &json,
        r#"{"a":{"chr1":"1-10,20-30"},"b":{"chr1":"5-25","chr2":"1-50"}}"#,
    )
    .unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["runlist", "combine", json.to_str().unwrap()])
        .run();
    assert_eq!(
        stdout,
        "{\n  \"chr1\": \"1-30\",\n  \"chr2\": \"1-50\"\n}\n"
    );
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "runlist",
            "combine",
            json.to_str().unwrap(),
            "--op",
            "intersect",
        ])
        .run();
    assert_eq!(
        stdout,
        "{\n  \"chr1\": \"5-10,20-25\",\n  \"chr2\": \"-\"\n}\n"
    );
}

#[test]
fn command_runlist_convert() {
    let dir = TempDir::new().unwrap();
    let json = dir.path().join("in.json");
    std::fs::write(&json, r#"{"chr1":"1-10,20-30"}"#).unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&["runlist", "convert", json.to_str().unwrap()])
        .run();
    assert_eq!(stdout, "chr1:1-10\nchr1:20-30\n");
    let (stdout, _) = PgrCmd::new()
        .args(&["runlist", "convert", json.to_str().unwrap(), "--longest"])
        .run();
    assert_eq!(stdout, "chr1:20-30\n");
}

#[test]
fn command_runlist_some() {
    let dir = TempDir::new().unwrap();
    let json = dir.path().join("in.json");
    let names = dir.path().join("names.txt");
    std::fs::write(&json, r#"{"chr1":"1-5","chr2":"6-9"}"#).unwrap();
    std::fs::write(&names, "chr1\n").unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "runlist",
            "some",
            json.to_str().unwrap(),
            names.to_str().unwrap(),
        ])
        .run();
    assert_eq!(stdout, "{\n  \"chr1\": \"1-5\"\n}\n");
}

#[test]
fn command_runlist_split() {
    let dir = TempDir::new().unwrap();
    let json = dir.path().join("multi.json");
    std::fs::write(&json, r#"{"a":{"chr1":"1-5"},"b":{"chr2":"6-9"}}"#).unwrap();
    let out = dir.path().join("out");
    PgrCmd::new()
        .args(&[
            "runlist",
            "split",
            json.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert_eq!(
        std::fs::read_to_string(out.join("a.json")).unwrap(),
        "{\"chr1\":\"1-5\"}\n"
    );
    assert_eq!(
        std::fs::read_to_string(out.join("b.json")).unwrap(),
        "{\"chr2\":\"6-9\"}\n"
    );
}

#[test]
fn command_runlist_stat() {
    let dir = TempDir::new().unwrap();
    let sizes = dir.path().join("sizes.txt");
    let json = dir.path().join("in.json");
    std::fs::write(&sizes, "chr1\t1000\nchr2\t2000\n").unwrap();
    std::fs::write(&json, r#"{"chr1":"1-500","chr2":"1-100"}"#).unwrap();
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "runlist",
            "stat",
            sizes.to_str().unwrap(),
            json.to_str().unwrap(),
        ])
        .run();
    assert_eq!(
        stdout,
        "chr,chrLength,size,coverage\n\
         chr1,1000,500,0.5000\n\
         chr2,2000,100,0.0500\n\
         all,3000,600,0.2000\n"
    );
}
