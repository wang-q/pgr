use assert_cmd::prelude::*;
use std::process::Command;

#[test]
fn command_indent() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("indent")
        .arg("tests/newick/hg38.7way.nwk")
        .arg("--text")
        .arg(".   ")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 19);
    assert!(stdout.contains(".   .   Human:"));
    assert!(stdout.contains("\n.   Opossum:"));

    Ok(())
}

#[test]
fn command_indent_compact() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("indent")
        .arg("tests/newick/catarrhini.nwk")
        .arg("--compact")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 1);
    assert_eq!(stdout.trim().lines().count(), 1); // Ensure only one line after trim
    assert!(stdout.contains("Gorilla"));

    Ok(())
}


