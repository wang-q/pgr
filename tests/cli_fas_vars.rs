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

#[test]
fn command_to_xlsx_many_sequences_no_overflow() -> anyhow::Result<()> {
    // A block with more than 32 ingroup sequences makes the per-variation
    // occurrence `pattern` a binary string wider than u32, which a direct
    // `from_str_radix` color-index parse used to fail on. It must succeed and
    // produce a valid workbook.
    let mut content = String::new();
    // First sequence differs at the last base to create a substitution.
    content.push_str(">A.chr1(+):1-4\nACGA\n");
    for i in 1..40 {
        content.push_str(&format!(">S{}.chr{}(+):1-4\nACGT\n", i, i));
    }
    let input = tempfile::NamedTempFile::new()?;
    std::io::Write::write_all(
        &mut std::fs::File::create(input.path())?,
        content.as_bytes(),
    )?;
    let out = tempfile::NamedTempFile::new()?;

    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("fas")
        .arg("to-xlsx")
        .arg(input.path())
        .arg("-o")
        .arg(out.path())
        .output()?;
    assert!(
        output.status.success(),
        "to-xlsx with >32 ingroup sequences must not fail: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let mut workbook: calamine::Xlsx<_> = calamine::open_workbook(out.path()).unwrap();
    let sheet = workbook.worksheet_range_at(0).unwrap().unwrap();
    // 40 ingroup rows for the single substitution plus a header row.
    assert!(
        sheet.height() >= 41,
        "unexpected sheet height: {}",
        sheet.height()
    );

    Ok(())
}

#[test]
fn command_to_xlsx_indel_fits_in_wrapped_section() -> anyhow::Result<()> {
    // A 3-base indel spans 3 columns. With --wrap 3 it fits exactly in one
    // section ending on the last (wrap) column. Regression for the
    // paint_indel off-by-one that needlessly wrapped such an indel into its
    // own section (and then wrapped again, leaving an empty section).
    let content = ">A.chr1(+):1-3\nAAA\n>B.chr1(+):1-3\n---\n";
    let input = tempfile::NamedTempFile::new()?;
    std::io::Write::write_all(
        &mut std::fs::File::create(input.path())?,
        content.as_bytes(),
    )?;
    let out = tempfile::NamedTempFile::new()?;

    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("fas")
        .arg("to-xlsx")
        .arg(input.path())
        .arg("--indel")
        .arg("--wrap")
        .arg("3")
        .arg("-o")
        .arg(out.path())
        .output()?;
    assert!(output.status.success());

    let mut wb: calamine::Xlsx<_> = calamine::open_workbook(out.path())?;
    let sheet = wb.worksheet_range_at(0).unwrap().unwrap();
    // sec_height = seq_count(2) + 1 + spacing(1) = 4. The single 3-base indel
    // fits entirely in section 1 (data at row 1, 0-based), so the cursor wraps
    // once and names appear in sections 1 and 2 (rows 1,2 and 5,6). There must
    // be no third, empty section (rows 9,10) — the off-by-one previously
    // wrapped the indel into its own section and then wrapped again, creating
    // one.
    assert!(
        sheet.get_value((1, 1)).is_some(),
        "indel should be drawn in section 1"
    );
    for row in [1u32, 2, 5, 6] {
        assert!(
            sheet.get_value((row, 0)).is_some(),
            "name should be present at row {}",
            row
        );
    }
    for row in [9u32, 10] {
        assert!(
            sheet.get_value((row, 0)).is_none(),
            "no empty third section: unexpected name at row {}",
            row
        );
    }

    Ok(())
}

#[test]
fn command_to_xlsx_colors_rejects_out_of_range() -> anyhow::Result<()> {
    // `--colors` outside [1, 15] must fail with a friendly error.
    let out = tempfile::NamedTempFile::new()?;
    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("fas")
        .arg("to-xlsx")
        .arg("tests/fas/example.fas")
        .arg("--colors")
        .arg("16")
        .arg("-o")
        .arg(out.path())
        .output()?;
    assert!(!output.status.success(), "--colors 16 must be rejected");
    let stderr = String::from_utf8(output.stderr)?;
    assert!(
        stderr.contains("--colors must be in [1, 15]"),
        "unexpected stderr: {}",
        stderr
    );

    Ok(())
}

#[test]
fn command_to_xlsx_colors_reduces_background_loop() -> anyhow::Result<()> {
    // With a small `--colors`, the per-variation background index is taken
    // modulo that count, so distinct variations reuse colors sooner. The
    // command must still succeed and produce a valid workbook.
    let out = tempfile::NamedTempFile::new()?;
    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("fas")
        .arg("to-xlsx")
        .arg("tests/fas/example.fas")
        .arg("--colors")
        .arg("3")
        .arg("-o")
        .arg(out.path())
        .output()?;
    assert!(
        output.status.success(),
        "to-xlsx --colors 3 must succeed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let mut workbook: calamine::Xlsx<_> = calamine::open_workbook(out.path()).unwrap();
    assert!(workbook.worksheet_range_at(0).is_some());

    Ok(())
}

#[test]
fn command_to_xlsx_spacing_and_wrapped_section_names() -> anyhow::Result<()> {
    // Two sequences differing at all 8 positions produce 8 substitutions.
    // With --wrap 3 each section holds 3 columns, so the variations span 3
    // sections. Every section must carry the sequence names, and --spacing
    // adds blank rows below each section (shifting later sections down).
    let mut content = String::new();
    content.push_str(">A.chr1(+):1-8\nAAAAAAAA\n>B.chr2(+):1-8\nCCCCCCCC\n");
    let input = tempfile::NamedTempFile::new()?;
    std::io::Write::write_all(
        &mut std::fs::File::create(input.path())?,
        content.as_bytes(),
    )?;

    let run = |spacing: &str, out: &std::path::Path| -> anyhow::Result<calamine::Xlsx<_>> {
        let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
        let output = cmd
            .arg("fas")
            .arg("to-xlsx")
            .arg(input.path())
            .arg("--wrap")
            .arg("3")
            .arg("--spacing")
            .arg(spacing)
            .arg("-o")
            .arg(out)
            .output()?;
        assert!(output.status.success());
        Ok(calamine::open_workbook(out).unwrap())
    };

    let out_default = tempfile::NamedTempFile::new()?;
    let mut wb_default = run("1", out_default.path())?;
    let default_sheet = wb_default.worksheet_range_at(0).unwrap().unwrap();

    let out_spaced = tempfile::NamedTempFile::new()?;
    let mut wb_spaced = run("3", out_spaced.path())?;
    let spaced_sheet = wb_spaced.worksheet_range_at(0).unwrap().unwrap();

    // section_height = seq_count(2) + 1 + spacing.
    // spacing=1 -> height 4: sections at rows 1-2, 5-6, 9-10 (names).
    // spacing=3 -> height 6: sections at rows 1-2, 7-8, 13-14 (names).
    for row in [2u32, 6, 10] {
        assert!(
            default_sheet.get_value((row, 0)).is_some(),
            "spacing=1 should name the section top at row {}",
            row
        );
    }
    for row in [2u32, 8, 14] {
        assert!(
            spaced_sheet.get_value((row, 0)).is_some(),
            "spacing=3 should name the section top at row {}",
            row
        );
    }

    Ok(())
}
