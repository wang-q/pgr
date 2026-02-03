use assert_cmd::prelude::*;
use std::process::Command;

#[test]
fn command_maf_to_fas() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("maf")
        .arg("to-fas")
        .arg("tests/maf/example.maf")
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();

    assert!(stdout.contains(">S288c.VIII(+):13377-13410"));
    assert!(stdout.contains("TTACTCGTCTTGCGGCCAAAACTCGAAGAAAAAC"));
    assert!(stdout.contains(">Spar.gi_29362578(-):72853-72885"));
    assert!(stdout.contains("TTACCCGTCTTGCGTCCAAAACTCGAA-AAAAAC"));
    assert_eq!(stdout.matches(">").count(), 8); // 2 blocks * 4 sequences
    assert_eq!(stdout.lines().count(), 18);
    assert!(stdout.contains("S288c.VIII"), "name list");
    assert!(stdout.contains(":42072-42168"), "coordinate transformed");

    Ok(())
}
