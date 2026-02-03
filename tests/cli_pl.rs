use assert_cmd::prelude::*;
use std::process::Command;
use tempfile::TempDir;

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
