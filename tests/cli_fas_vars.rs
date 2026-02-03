use assert_cmd::prelude::*;
use calamine::Reader;
use std::process::Command;
use tempfile::NamedTempFile;


#[test]
fn command_variation() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fas")
        .arg("variation")
        .arg("tests/fas/example.fas")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 81);

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fas")
        .arg("variation")
        .arg("tests/fas/example.fas")
        .arg("--outgroup")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 49);

    Ok(())
}

#[test]
fn command_toxlsx() -> anyhow::Result<()> {
    let temp_file = NamedTempFile::new()?.into_temp_path();
    let temp_path = temp_file.to_str().unwrap();

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fas")
        .arg("toxlsx")
        .arg("tests/fas/example.fas")
        .arg("-o")
        .arg(temp_path)
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 0);
    let mut workbook: calamine::Xlsx<_> = calamine::open_workbook(temp_path).unwrap();
    let sheet = workbook.worksheet_range_at(0).unwrap().unwrap();

    // row-col
    assert_eq!(
        sheet.get_value((1, 1)).unwrap().to_string(),
        "G".to_string()
    );
    assert_eq!(
        sheet.get_value((19, 8)).unwrap().to_string(),
        "C".to_string()
    );

    Ok(())
}

#[test]
fn command_toxlsx_indel() -> anyhow::Result<()> {
    let temp_file = NamedTempFile::new()?.into_temp_path();
    let temp_path = temp_file.to_str().unwrap();

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fas")
        .arg("toxlsx")
        .arg("tests/fas/example.fas")
        .arg("--indel")
        .arg("-o")
        .arg(temp_path)
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;
    assert_eq!(stdout.lines().count(), 0);

    let mut workbook: calamine::Xlsx<_> = calamine::open_workbook(temp_path).unwrap();
    let sheet = workbook.worksheet_range_at(0).unwrap().unwrap();

    assert_eq!(
        sheet.get_value((1, 1)).unwrap().to_string(),
        "G".to_string()
    );
    assert_eq!(
        sheet.get_value((19, 8)).unwrap().to_string(),
        "D1".to_string()
    );

    Ok(())
}

#[test]
fn command_toxlsx_nocomplex() -> anyhow::Result<()> {
    let temp_file = NamedTempFile::new()?.into_temp_path();
    let temp_path = temp_file.to_str().unwrap();

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fas")
        .arg("toxlsx")
        .arg("tests/fas/example.fas")
        .arg("--indel")
        .arg("--nocomplex")
        .arg("-o")
        .arg(temp_path)
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;
    assert_eq!(stdout.lines().count(), 0);

    let mut workbook: calamine::Xlsx<_> = calamine::open_workbook(temp_path).unwrap();
    let sheet = workbook.worksheet_range_at(0).unwrap().unwrap();

    assert_eq!(
        sheet.get_value((13, 7)).unwrap().to_string(),
        "D1".to_string()
    );
    assert_eq!(
        sheet.get_value((13, 8)).unwrap().to_string(),
        "T".to_string()
    );

    Ok(())
}

#[test]
fn command_toxlsx_nosingle() -> anyhow::Result<()> {
    let temp_file = NamedTempFile::new()?.into_temp_path();
    let temp_path = temp_file.to_str().unwrap();

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fas")
        .arg("toxlsx")
        .arg("tests/fas/example.fas")
        .arg("--indel")
        .arg("--nosingle")
        .arg("-o")
        .arg(temp_path)
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;
    assert_eq!(stdout.lines().count(), 0);

    let mut workbook: calamine::Xlsx<_> = calamine::open_workbook(temp_path).unwrap();
    let sheet = workbook.worksheet_range_at(0).unwrap().unwrap();

    assert_eq!(
        sheet.get_value((13, 3)).unwrap().to_string(),
        "I1".to_string()
    );
    assert_eq!(
        sheet.get_value((13, 4)).unwrap().to_string(),
        "G".to_string()
    );

    Ok(())
}

#[test]
fn command_toxlsx_minmax() -> anyhow::Result<()> {
    let temp_file = NamedTempFile::new()?.into_temp_path();
    let temp_path = temp_file.to_str().unwrap();

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fas")
        .arg("toxlsx")
        .arg("tests/fas/example.fas")
        .arg("--indel")
        .arg("--min")
        .arg("0.3")
        .arg("--max")
        .arg("0.7")
        .arg("-o")
        .arg(temp_path)
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;
    assert_eq!(stdout.lines().count(), 0);

    let mut workbook: calamine::Xlsx<_> = calamine::open_workbook(temp_path).unwrap();
    let sheet = workbook.worksheet_range_at(0).unwrap().unwrap();

    assert_eq!(
        sheet.get_value((13, 1)).unwrap().to_string(),
        "D1".to_string()
    );
    assert_eq!(
        sheet.get_value((13, 5)).unwrap().to_string(),
        "T".to_string()
    );

    Ok(())
}

#[test]
fn command_toxlsx_outgroup() -> anyhow::Result<()> {
    let temp_file = NamedTempFile::new()?.into_temp_path();
    let temp_path = temp_file.to_str().unwrap();

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fas")
        .arg("toxlsx")
        .arg("tests/fas/example.fas")
        .arg("--indel")
        .arg("--outgroup")
        .arg("-o")
        .arg(temp_path)
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;
    assert_eq!(stdout.lines().count(), 0);

    let mut workbook: calamine::Xlsx<_> = calamine::open_workbook(temp_path).unwrap();
    let sheet = workbook.worksheet_range_at(0).unwrap().unwrap();

    assert_eq!(
        sheet.get_value((7, 1)).unwrap().to_string(),
        "A".to_string()
    );
    assert_eq!(
        sheet.get_value((14, 4)).unwrap().to_string(),
        "I1".to_string()
    );

    Ok(())
}
