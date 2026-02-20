use assert_cmd::cargo::cargo_bin_cmd;
use std::fs;
use tempfile::TempDir;

#[test]
fn command_fas_multiz_core() -> anyhow::Result<()> {
    let tempdir = TempDir::new()?;
    let out_path = tempdir.path().join("merged.fas");
    let out_str = out_path.to_str().unwrap();

    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("fas")
        .arg("multiz")
        .arg("-r")
        .arg("S288c")
        .arg("tests/fas/S288cvsRM11_1a.slice.fas")
        .arg("tests/fas/S288cvsYJM789.slice.fas")
        .arg("tests/fas/S288cvsSpar.slice.fas")
        .arg("-o")
        .arg(out_str)
        .output()?;

    if !output.status.success() {
        let stderr = String::from_utf8(output.stderr)?;
        panic!("pgr fas multiz failed: {}", stderr);
    }

    assert!(out_path.is_file());
    let content = fs::read_to_string(out_path)?;
    assert!(content.lines().count() > 0);

    tempdir.close()?;

    Ok(())
}

#[test]
fn command_fas_multiz_affine_gap() -> anyhow::Result<()> {
    let tempdir = TempDir::new()?;
    let out_path = tempdir.path().join("merged_affine.fas");
    let out_str = out_path.to_str().unwrap();

    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("fas")
        .arg("multiz")
        .arg("-r")
        .arg("S288c")
        .arg("tests/fas/S288cvsRM11_1a.slice.fas")
        .arg("tests/fas/S288cvsYJM789.slice.fas")
        .arg("tests/fas/S288cvsSpar.slice.fas")
        .arg("--gap-open")
        .arg("400")
        .arg("--gap-extend")
        .arg("30")
        .arg("-o")
        .arg(out_str)
        .output()?;

    if !output.status.success() {
        let stderr = String::from_utf8(output.stderr)?;
        panic!("pgr fas multiz (affine) failed: {}", stderr);
    }

    assert!(out_path.is_file());
    let content = fs::read_to_string(out_path)?;
    assert!(content.lines().count() > 0);

    tempdir.close()?;

    Ok(())
}

#[test]
fn command_fas_multiz_custom_matrix() -> anyhow::Result<()> {
    let tempdir = TempDir::new()?;
    let out_path = tempdir.path().join("merged_matrix.fas");
    let out_str = out_path.to_str().unwrap();

    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("fas")
        .arg("multiz")
        .arg("-r")
        .arg("S288c")
        .arg("tests/fas/S288cvsRM11_1a.slice.fas")
        .arg("tests/fas/S288cvsYJM789.slice.fas")
        .arg("tests/fas/S288cvsSpar.slice.fas")
        .arg("--score-matrix")
        .arg("hoxd55")
        .arg("-o")
        .arg(out_str)
        .output()?;

    if !output.status.success() {
        let stderr = String::from_utf8(output.stderr)?;
        panic!("pgr fas multiz (custom matrix) failed: {}", stderr);
    }

    assert!(out_path.is_file());
    let content = fs::read_to_string(out_path)?;
    assert!(content.lines().count() > 0);

    tempdir.close()?;

    Ok(())
}
