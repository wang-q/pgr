use assert_cmd::prelude::*;
use std::process::Command;
use tempfile::TempDir;

#[test]
fn command_pl_prefilter_help() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("pl").arg("prefilter").arg("--help").output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("Prefilter genome/metagenome assembly by amino acid minimizers"));
    Ok(())
}

#[test]
fn command_pl_prefilter_run() -> anyhow::Result<()> {
    let tempdir = TempDir::new()?;
    let _tempdir_str = tempdir.path().to_str().unwrap();

    let input = "tests/fas/NC_000932.fa";
    
    // Generate protein reference
    let ref_file = tempdir.path().join("ref.pep.fa");
    let ref_path = ref_file.to_str().unwrap();
    
    let mut cmd_gen = Command::cargo_bin("pgr")?;
    let output_gen = cmd_gen
        .arg("fa")
        .arg("six-frame")
        .arg(input)
        .arg("--outfile")
        .arg(ref_path)
        .output()?;
    assert!(output_gen.status.success());

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("pl")
        .arg("prefilter")
        .arg(input)
        .arg(ref_path)
        .arg("--chunk")
        .arg("50000")
        .output()?;
    
    // Check for success
    if !output.status.success() {
        let stderr = String::from_utf8(output.stderr)?;
        println!("stderr: {}", stderr);
    }
    assert!(output.status.success());

    let stdout = String::from_utf8(output.stdout)?;
    // Since we compare the file against its translation, we should get matches.
    assert!(!stdout.is_empty());

    Ok(())
}

#[test]
fn command_pl_p2m() -> anyhow::Result<()> {
    match which::which("spanr") {
        Err(_) => return Ok(()),
        Ok(_) => {}
    }

    let tempdir = TempDir::new()?;
    let tempdir_str = tempdir.path().to_str().unwrap();

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("pl")
        .arg("p2m")
        .arg("tests/fas/S288cvsRM11_1a.slice.fas")
        .arg("tests/fas/S288cvsYJM789.slice.fas")
        .arg("tests/fas/S288cvsSpar.slice.fas")
        .arg("-o")
        .arg(tempdir_str)
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 13);
    assert!(&tempdir.path().join("merge.json").is_file());
    assert!(&tempdir.path().join("join.subset.fas").is_file());

    tempdir.close()?;

    Ok(())
}

#[test]
fn command_pl_trf_help() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("pl").arg("trf").arg("--help").output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("Identify tandem repeats in a genome"));
    Ok(())
}

#[test]
fn command_pl_ir_help() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("pl").arg("ir").arg("--help").output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("Identify interspersed repeats in a genome"));
    Ok(())
}

#[test]
fn command_pl_rept_help() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("pl").arg("rept").arg("--help").output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("Identify repetitive regions in a genome"));
    Ok(())
}

#[test]
fn command_pl_ucsc_help() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("pl").arg("ucsc").arg("--help").output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("UCSC chain/net pipeline"));
    Ok(())
}
