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
fn command_refine_msa() -> anyhow::Result<()> {
    let mut bin = String::new();
    for e in &["clustalw", "clustal-w", "clustalw2"] {
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
        .arg("refine")
        .arg("tests/fas/refine.fas")
        .arg("--msa")
        .arg("clustalw")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 18);
    assert!(stdout.contains("---"), "dashes added");

    // --outgroup
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fas")
        .arg("refine")
        .arg("tests/fas/refine2.fas")
        .arg("--msa")
        .arg("clustalw")
        .arg("--outgroup")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 7);
    assert!(stdout.contains("CA-GT"), "outgroup trimmed");

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fas")
        .arg("refine")
        .arg("tests/fas/refine2.fas")
        .arg("--msa")
        .arg("clustalw")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 7);
    assert!(stdout.contains("CA--GT"), "outgroup not trimmed");

    // quick
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fas")
        .arg("refine")
        .arg("tests/fas/refine2.fas")
        .arg("--msa")
        .arg("clustalw")
        .arg("--quick")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 7);
    assert!(stdout.contains("CA--GT"), "outgroup not trimmed");

    Ok(())
}

#[test]
fn command_refine_poa() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fas")
        .arg("refine")
        .arg("tests/fas/refine.fas")
        .arg("--msa")
        .arg("builtin")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 18);
    // The exact alignment might vary, but it should contain dashes
    assert!(stdout.contains("-"), "dashes added by builtin");
    // Check if sequences are still present
    assert!(stdout.contains(">S288c"));
    assert!(stdout.contains(">Spar"));

    Ok(())
}

#[test]
fn command_refine_default() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fas")
        .arg("refine")
        .arg("tests/fas/refine.fas")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 18);
    // Should be builtin by default, so it should contain dashes
    assert!(stdout.contains("-"), "dashes added by default builtin");
    assert!(stdout.contains(">S288c"));
    assert!(stdout.contains(">Spar"));

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

#[test]
fn command_consensus_comparison() -> anyhow::Result<()> {
    // Check if spoa is available
    if which::which("spoa").is_err() {
        return Ok(());
    }

    let mut cmd_builtin = Command::cargo_bin("pgr")?;
    let output_builtin = cmd_builtin
        .arg("fas")
        .arg("consensus")
        .arg("tests/fas/refine.fas")
        .arg("--engine")
        .arg("builtin")
        .output()?;
    let stdout_builtin = String::from_utf8(output_builtin.stdout)?;

    let mut cmd_spoa = Command::cargo_bin("pgr")?;
    let output_spoa = cmd_spoa
        .arg("fas")
        .arg("consensus")
        .arg("tests/fas/refine.fas")
        .arg("--engine")
        .arg("spoa")
        .output()?;
    let stdout_spoa = String::from_utf8(output_spoa.stdout)?;

    assert_eq!(stdout_builtin, stdout_spoa, "Builtin and Spoa outputs should match");

    Ok(())
}

#[test]
fn command_consensus_options() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fas")
        .arg("consensus")
        .arg("tests/fas/refine.fas")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 6);
    assert!(stdout.contains(">consensus\n"), "simple name");
    assert!(stdout.contains(">consensus.I("), "fas name");

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fas")
        .arg("consensus")
        .arg("tests/fas/refine.fas")
        .arg("--outgroup")
        .arg("--parallel")
        .arg("2")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 10);
    assert!(stdout.contains(">Spar"), "outgroup");

    Ok(())
}
