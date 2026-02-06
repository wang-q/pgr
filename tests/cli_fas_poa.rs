use assert_cmd::prelude::*;
use std::process::Command;

#[test]
fn command_consensus_builtin() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fas")
        .arg("consensus")
        .arg("tests/fas/refine.fas")
        .arg("--engine")
        .arg("builtin")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 6);
    assert!(stdout.contains(">consensus\n"), "simple name");
    assert!(stdout.contains(">consensus.I("), "fas name");

    Ok(())
}

#[test]
fn command_consensus_spoa() -> anyhow::Result<()> {
    let mut bin = String::new();
    for e in &["spoa"] {
        if let Ok(pth) = which::which(e) {
            bin = pth.to_string_lossy().to_string();
            break;
        }
    }
    if bin.is_empty() {
        return Ok(());
    } else {
        eprintln!("bin = {:#?}", bin);
    }

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fas")
        .arg("consensus")
        .arg("tests/fas/refine.fas")
        .arg("--engine")
        .arg("spoa")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 6);
    assert!(stdout.contains(">consensus\n"), "simple name");
    assert!(stdout.contains(">consensus.I("), "fas name");

    Ok(())
}

#[test]
fn command_consensus_params() -> anyhow::Result<()> {
    // Test with custom parameters
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fas")
        .arg("consensus")
        .arg("tests/fas/refine.fas")
        .arg("--match")
        .arg("2")
        .arg("--mismatch")
        .arg("-3")
        .arg("--gap-open")
        .arg("-5")
        .arg("--gap-extend")
        .arg("-1")
        .arg("--algorithm")
        .arg("global")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;
    let stderr = String::from_utf8(output.stderr)?;

    if !output.status.success() {
        println!("stderr: {}", stderr);
    }

    assert_eq!(stdout.lines().count(), 6);
    assert!(stdout.contains(">consensus\n"), "simple name");

    Ok(())
}
