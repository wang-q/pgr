use assert_cmd::prelude::*;
use std::process::Command;

#[test]
fn command_stat() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("stat")
        .arg("tests/newick/hg38.7way.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 5);
    assert!(stdout.contains("leaf labels\t7"));

    Ok(())
}
