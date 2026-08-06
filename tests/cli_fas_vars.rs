use calamine::Reader;
use tempfile::NamedTempFile;

#[test]
fn command_variation() -> anyhow::Result<()> {
    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("fas")
        .arg("variation")
        .arg("tests/fas/example.fas")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 81);

    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
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
fn command_to_xlsx() -> anyhow::Result<()> {
    let temp_file = NamedTempFile::new()?.into_temp_path();
    let temp_path = temp_file.to_str().unwrap();

    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("fas")
        .arg("to-xlsx")
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
fn command_to_xlsx_indel() -> anyhow::Result<()> {
    let temp_file = NamedTempFile::new()?.into_temp_path();
    let temp_path = temp_file.to_str().unwrap();

    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("fas")
        .arg("to-xlsx")
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
fn command_to_xlsx_nocomplex() -> anyhow::Result<()> {
    let temp_file = NamedTempFile::new()?.into_temp_path();
    let temp_path = temp_file.to_str().unwrap();

    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("fas")
        .arg("to-xlsx")
        .arg("tests/fas/example.fas")
        .arg("--indel")
        .arg("--no-complex")
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
fn command_to_xlsx_nosingle() -> anyhow::Result<()> {
    let temp_file = NamedTempFile::new()?.into_temp_path();
    let temp_path = temp_file.to_str().unwrap();

    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("fas")
        .arg("to-xlsx")
        .arg("tests/fas/example.fas")
        .arg("--indel")
        .arg("--no-single")
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
fn command_to_xlsx_minmax() -> anyhow::Result<()> {
    let temp_file = NamedTempFile::new()?.into_temp_path();
    let temp_path = temp_file.to_str().unwrap();

    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("fas")
        .arg("to-xlsx")
        .arg("tests/fas/example.fas")
        .arg("--indel")
        .arg("--min-freq")
        .arg("0.3")
        .arg("--max-freq")
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
fn command_to_xlsx_outgroup() -> anyhow::Result<()> {
    let temp_file = NamedTempFile::new()?.into_temp_path();
    let temp_path = temp_file.to_str().unwrap();

    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("fas")
        .arg("to-xlsx")
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

#[test]
fn command_to_xlsx_outgroup_ambiguity_no_error() -> anyhow::Result<()> {
    // The outgroup (last sequence) has an IUPAC ambiguity code ('R') at a
    // position where the ingroup sequences are polymorphic canonical bases.
    // The outgroup-cell color style for 'R' is not registered, so the export
    // must fall back to a default style instead of failing.
    let input = tempfile::NamedTempFile::new()?;
    std::io::Write::write_all(
        &mut std::fs::File::create(input.path())?,
        b">A.chr1(+):1-10\nAAAATTTTGG\n>B.chr2(+):1-10\nAAAATTTTAG\n>C.chr3(+):1-10\nAAAATTTTRG\n",
    )?;
    let out = tempfile::NamedTempFile::new()?;

    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("fas")
        .arg("to-xlsx")
        .arg(input.path())
        .arg("--outgroup")
        .arg("-o")
        .arg(out.path())
        .output()?;
    assert!(
        output.status.success(),
        "to-xlsx must succeed with an ambiguous outgroup: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    Ok(())
}

#[test]
fn command_to_xlsx_shared_all_gap_region_no_error() -> anyhow::Result<()> {
    // A block where every sequence shares an all-gap region (here positions
    // 2-3) has no variation in that span. `to-xlsx --indel` must skip such a
    // span instead of bailing and aborting the whole command.
    let input = tempfile::NamedTempFile::new()?;
    std::io::Write::write_all(
        &mut std::fs::File::create(input.path())?,
        b">A.chr1(+):1-4\nA--A\n>B.chr2(+):1-4\nA--A\n",
    )?;
    let out = tempfile::NamedTempFile::new()?;

    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("fas")
        .arg("to-xlsx")
        .arg(input.path())
        .arg("--indel")
        .arg("-o")
        .arg(out.path())
        .output()?;
    assert!(
        output.status.success(),
        "to-xlsx must skip a shared all-gap region: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    Ok(())
}

#[test]
fn command_variation_outgroup_unequal_length_no_panic() -> anyhow::Result<()> {
    // Last (outgroup) sequence is shorter than the ingroup sequences: the
    // polarization step must return a friendly error instead of panicking.
    let malformed = tempfile::NamedTempFile::new()?;
    std::io::Write::write_all(
        &mut std::fs::File::create(malformed.path())?,
        b">A.chr1(+):1-10\nAAAATTTTGG\n>B.chr2(+):1-10\nAAAATTTTAG\n>C.chr3(+):1-8\nAAAATTTT\n",
    )?;

    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("fas")
        .arg("variation")
        .arg(malformed.path())
        .arg("--outgroup")
        .output()?;
    assert!(
        !output.status.success(),
        "expected a friendly error, not a panic"
    );
    let stderr = String::from_utf8(output.stderr)?;
    assert!(
        stderr.contains("outgroup sequence too short"),
        "unexpected stderr: {}",
        stderr
    );

    Ok(())
}

#[test]
fn command_to_xlsx_unequal_length_no_panic() -> anyhow::Result<()> {
    // Varying sequence lengths within a block must error gracefully, not panic.
    let malformed = tempfile::NamedTempFile::new()?;
    std::io::Write::write_all(
        &mut std::fs::File::create(malformed.path())?,
        b">A.chr1(+):1-10\nAAAATTTTGG\n>B.chr2(+):1-10\nAAAATTTTAG\n>C.chr3(+):1-8\nAAAATTTT\n",
    )?;
    let out = tempfile::NamedTempFile::new()?;

    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("fas")
        .arg("to-xlsx")
        .arg(malformed.path())
        .arg("--indel")
        .arg("--outgroup")
        .arg("-o")
        .arg(out.path())
        .output()?;
    assert!(
        !output.status.success(),
        "expected a friendly error, not a panic"
    );
    let stderr = String::from_utf8(output.stderr)?;
    assert!(
        stderr.contains("unequal lengths") || stderr.contains("outgroup sequence too short"),
        "unexpected stderr: {}",
        stderr
    );

    Ok(())
}
