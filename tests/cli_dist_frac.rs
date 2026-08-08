#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::path::PathBuf;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/dist/input")
        .join(name)
}

#[test]
fn command_dist_frac_basic() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "frac",
            fixture("random.fa").to_str().unwrap(),
            "-k",
            "21",
            "--scale",
            "10",
            "--zero",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 9);
    assert!(stdout.contains("r1\tr1\t0.0000\t1.0000\t1.0000"));
    assert!(stdout.contains("r1\tr2\t0.1365\t0.0293\t0.0531"));
    assert!(stdout.contains("r1\tr3\t1.0000\t0.0000\t0.0000"));
}

#[test]
fn command_dist_frac_ci() {
    // --ci appends the 95% ANI confidence interval (2 extra fields).
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "frac",
            fixture("random.fa").to_str().unwrap(),
            "-k",
            "21",
            "--scale",
            "10",
            "--zero",
            "--ci",
        ])
        .run();

    let line = stdout.lines().find(|l| l.starts_with("r1\tr2\t")).unwrap();
    let fields: Vec<&str> = line.split_whitespace().collect();
    assert_eq!(fields.len(), 7, "ci output fields: {line}");
    // ANI = 1 - mash = 1 - 0.1365 = 0.8635; CI brackets it.
    let ani: f64 = 1.0 - fields[2].parse::<f64>().unwrap();
    let lo: f64 = fields[5].parse().unwrap();
    let hi: f64 = fields[6].parse().unwrap();
    assert!(lo < ani && hi > ani, "CI must bracket ANI: {line}");
}

#[test]
fn command_dist_frac_sim() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "frac",
            fixture("random.fa").to_str().unwrap(),
            "-k",
            "21",
            "--scale",
            "10",
            "--zero",
            "--sim",
        ])
        .run();

    assert!(stdout.contains("r1\tr1\t1.0000\t1.0000\t1.0000"));
}
