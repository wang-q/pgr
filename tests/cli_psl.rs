use assert_cmd::Command;
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

fn get_path(subcommand: &str, dir: &str, filename: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("tests/psl");
    path.push(subcommand);
    path.push(dir);
    path.push(filename);
    path
}

//
// psl histo
//

#[test]
fn test_histo_apq_base() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("histo", "input", "basic.psl");
    let output = temp.path().join("apq.histo");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("psl")
        .arg("histo")
        .arg("--what")
        .arg("alignsPerQuery")
        .arg(&input)
        .arg("-o")
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
fn test_histo_apq_multi() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("histo", "input", "basic.psl");
    let output = temp.path().join("apq_multi.histo");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("psl")
        .arg("histo")
        .arg("--what")
        .arg("alignsPerQuery")
        .arg(&input)
        .arg("-o")
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
fn test_histo_cover_spread() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("histo", "input", "basic.psl");
    let output = temp.path().join("cover.histo");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("psl")
        .arg("histo")
        .arg("--what")
        .arg("coverSpread")
        .arg(&input)
        .arg("-o")
        .arg(&output);
    cmd.assert().success();

    // NM_000014.3: 1 align. Spread = 0.
    // NM_005577.1: 4 aligns.
    //   3335+96+0 / 13938 = 0.24616
    //   3444+105+0 / 13938 = 0.25463
    //   3482+120+0 / 13938 = 0.25843
    //   6410+4+0 / 13938 = 0.46018
    //   Diff: 0.46018 - 0.24616 = 0.2140
    //
    // Just checking it runs and produces output. Precise float matching is tricky.
    // I will check if output contains "0.2140"
    let output_content = fs::read_to_string(&output)?;
    assert!(output_content.contains("0.2140"));
    assert!(output_content.contains("0.0000")); // Singletons or identicals

    Ok(())
}

//
// psl tochain
//

#[test]
fn test_tochain_fix_strand() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("tochain", "input", "mtor.psl");
    let expected_output = get_path("tochain", "expected", "example3.chain");
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
fn test_tochain_fail_neg_strand() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("tochain", "input", "mtor.psl");
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

//
// psl rc
//

#[test]
fn test_rc_mrna() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("psl")
        .arg("rc")
        .arg(get_path("rc", "input", "mrna.psl"))
        .arg("-o")
        .arg("stdout")
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();
    let expected = std::fs::read_to_string(get_path("rc", "expected", "mrnaTest.psl"))?;

    assert_eq!(stdout.replace("\r\n", "\n"), expected.replace("\r\n", "\n"));
    Ok(())
}

#[test]
fn test_rc_trans() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("psl")
        .arg("rc")
        .arg(get_path("rc", "input", "trans.psl"))
        .arg("-o")
        .arg("stdout")
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();
    let expected = std::fs::read_to_string(get_path("rc", "expected", "transTest.psl"))?;

    assert_eq!(stdout.replace("\r\n", "\n"), expected.replace("\r\n", "\n"));
    Ok(())
}

//
// psl swap
//

#[test]
fn test_swap_mrna() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("psl")
        .arg("swap")
        .arg(get_path("swap", "input", "mrna.psl"))
        .arg("-o")
        .arg("stdout")
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();
    let expected = std::fs::read_to_string(get_path("swap", "expected", "mrnaTest.psl"))?;

    assert_eq!(stdout.replace("\r\n", "\n"), expected.replace("\r\n", "\n"));
    Ok(())
}

#[test]
fn test_swap_trans() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("psl")
        .arg("swap")
        .arg(get_path("swap", "input", "trans.psl"))
        .arg("-o")
        .arg("stdout")
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();
    let expected = std::fs::read_to_string(get_path("swap", "expected", "transTest.psl"))?;

    assert_eq!(stdout.replace("\r\n", "\n"), expected.replace("\r\n", "\n"));
    Ok(())
}

#[test]
fn test_swap_mrna_no_rc() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("psl")
        .arg("swap")
        .arg("--no-rc")
        .arg(get_path("swap", "input", "mrna.psl"))
        .arg("-o")
        .arg("stdout")
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();
    let expected = std::fs::read_to_string(get_path("swap", "expected", "mrnaNoRcTest.psl"))?;

    assert_eq!(stdout.replace("\r\n", "\n"), expected.replace("\r\n", "\n"));
    Ok(())
}
