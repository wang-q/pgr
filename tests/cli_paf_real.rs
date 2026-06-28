#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

// ── real-data validation ────────────────────────────────────────

#[test]
fn command_maf_to_paf_real_multiz_spar() {
    let (stdout, _) = PgrCmd::new()
        .args(&["maf", "to-paf", "tests/multiz/S288cvsSpar.maf"])
        .run();
    let fields: Vec<&str> = stdout.trim().split('\t').collect();
    assert_eq!(fields[0], "Spar.gi_29362594");
    assert_eq!(fields[4], "-");
    assert_eq!(fields[5], "S288c.I");
}

#[test]
fn command_maf_to_paf_real_multiz_rm11() {
    let (stdout, _) = PgrCmd::new()
        .args(&["maf", "to-paf", "tests/multiz/S288cvsRM11_1a.maf"])
        .run();
    let fields: Vec<&str> = stdout.trim().split('\t').collect();
    assert_eq!(fields[0], "RM11_1a.scaffold_17");
    assert_eq!(fields[4], "+");
    assert_eq!(fields[5], "S288c.I");
}

#[test]
fn command_paf_query_real_multiz_transitive() {
    use std::fs;
    use std::process::Command;
    let temp = tempfile::TempDir::new().unwrap();
    let paf_path = temp.path().join("merged.paf");
    let idx_path = temp.path().join("merged.paf.idx");
    let pgr = std::env::current_dir().unwrap().join("target/debug/pgr");
    let spar_out = Command::new(&pgr)
        .args(["maf", "to-paf", "tests/multiz/S288cvsSpar.maf"])
        .output()
        .unwrap();
    let rm11_out = Command::new(&pgr)
        .args(["maf", "to-paf", "tests/multiz/S288cvsRM11_1a.maf"])
        .output()
        .unwrap();
    let mut merged = spar_out.stdout.clone();
    merged.extend_from_slice(&rm11_out.stdout);
    fs::write(&paf_path, &merged).unwrap();
    let (_, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "index",
            paf_path.to_str().unwrap(),
            "-o",
            idx_path.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("saved to"));
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "paf",
            "query",
            idx_path.to_str().unwrap(),
            "S288c.I:26000-30000",
            "--transitive",
        ])
        .run();
    assert!(stderr.contains("Loading index"));
    assert!(stdout.contains("Spar.gi_29362594"), "Spar not found");
    assert!(stdout.contains("RM11_1a.scaffold_17"), "RM11 not found");
}
