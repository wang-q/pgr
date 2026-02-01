use assert_cmd::Command;
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

#[test]
fn command_psl_histo_apq_base() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/psl/basic.psl");
    let output = temp.path().join("apq.histo");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("psl")
        .arg("histo")
        .arg("alignsPerQuery")
        .arg(&input)
        .arg(&output);
    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    // Check for expected counts.
    // NM_033178.1: 2
    // NM_173571.1: 3
    // NM_000014.3: 1
    // NM_000015.1: 1
    // NM_153248.2: 1
    // NM_005577.1: 4
    // NM_FAKE.1: 2
    // Expected output order depends on hash map iteration unless sorted. 
    // I implemented sorting by key.
    // Sorted keys: NM_000014.3, NM_000015.1, NM_005577.1, NM_033178.1, NM_153248.2, NM_173571.1, NM_FAKE.1
    // Counts: 1, 1, 4, 2, 1, 3, 2
    let expected = "1\n1\n4\n2\n1\n3\n2\n";
    assert_eq!(output_content, expected);

    Ok(())
}

#[test]
fn command_psl_histo_apq_multi() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/psl/basic.psl");
    let output = temp.path().join("apq_multi.histo");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("psl")
        .arg("histo")
        .arg("alignsPerQuery")
        .arg(&input)
        .arg(&output)
        .arg("--multi-only");
    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    // Multi only: NM_005577.1 (4), NM_033178.1 (2), NM_173571.1 (3), NM_FAKE.1 (2)
    // Order: NM_005577.1, NM_033178.1, NM_173571.1, NM_FAKE.1
    // Counts: 4, 2, 3, 2
    let expected = "4\n2\n3\n2\n";
    assert_eq!(output_content, expected);

    Ok(())
}

#[test]
fn command_psl_histo_cover_spread() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/psl/basic.psl");
    let output = temp.path().join("cover.histo");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("psl")
        .arg("histo")
        .arg("coverSpread")
        .arg(&input)
        .arg(&output);
    cmd.assert().success();

    // NM_000014.3: 1 align. Spread = 0.
    // NM_005577.1: 4 aligns. 
    //   3335+96+0 / 13938 = 0.24616
    //   3444+105+0 / 13938 = 0.25463
    //   3482+120+0 / 13938 = 0.25843
    //   6410+4+0 / 13938 = 0.46018
    //   Diff: 0.46018 - 0.24616 = 0.2140
    
    // Just checking it runs and produces output. Precise float matching is tricky.
    // I will check if output contains "0.2140"
    let output_content = fs::read_to_string(&output)?;
    assert!(output_content.contains("0.2140"));
    assert!(output_content.contains("0.0000")); // Singletons or identicals

    Ok(())
}

#[test]
fn command_psl_to_chain_fix_strand() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/psl/mtor.psl");
    let expected_output = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/psl/example3.chain");
    let output = temp.path().join("output.chain");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("psl")
        .arg("tochain")
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .arg("--fix-strand");
    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    let expected_content = fs::read_to_string(&expected_output)?;
    
    assert_eq!(output_content, expected_content);

    Ok(())
}

#[test]
fn command_psl_to_chain_fail_neg_strand() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/psl/mtor.psl");
    let output = temp.path().join("fail.chain");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("psl")
        .arg("tochain")
        .arg(&input)
        .arg("-o")
        .arg(&output);
    // Should fail because mtor.psl has negative target strand
    cmd.assert().failure();

    Ok(())
}
