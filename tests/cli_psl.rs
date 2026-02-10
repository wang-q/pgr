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
// psl to-chain
//

#[test]
fn test_to_chain_fix_strand() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("to_chain", "input", "mtor.psl");
    let expected_output = get_path("to_chain", "expected", "example3.chain");
    let output = temp.path().join("out.chain");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("psl")
        .arg("to-chain")
        .arg(&input)
        .arg("--output")
        .arg(&output)
        .arg("--fix-strand");
    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    let expected_content = fs::read_to_string(&expected_output)?;
    assert_eq!(output_content, expected_content);

    Ok(())
}

#[test]
fn test_to_chain_fail_neg_strand() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("to_chain", "input", "mtor.psl");
    let output = temp.path().join("out.chain");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("psl")
        .arg("to-chain")
        .arg(&input)
        .arg("-o")
        .arg(&output);

    // Should fail because of negative target strand without --fix-strand
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
// psl lift
//

#[test]
fn test_lift_basic() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("lift", "", "test_fragment.psl");
    let sizes = get_path("lift", "", "chrom.sizes");
    let output = temp.path().join("lifted.psl");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("psl")
        .arg("lift")
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .arg("--q-sizes")
        .arg(&sizes);
    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    
    // Expected output check
    // The input file contains two records.
    // First record: chr1:101-200 (+), qStart=10, qEnd=20 (on fragment).
    //   Lifted: chr1 (+), qStart=100+10=110, qEnd=100+20=120.
    // Second record: chr1:101-200 (-), qStart=10, qEnd=20 (on RC fragment).
    //   Lifted: chr1 (-), qStart=810, qEnd=820.
    //   Calculation for (-):
    //     q_size=1000. offset=100.
    //     new_qStart = 1000 - (100 - 10 + 100) = 810.
    //     new_qEnd = 1000 - (100 - 20 + 100) = 820.
    
    // Check first record
    assert!(output_content.contains("chr1\t1000\t110\t120"));
    // Check second record
    assert!(output_content.contains("chr1\t1000\t810\t820"));
    // Check that qStarts for blocks are also correct
    // First record block start: 110
    assert!(output_content.contains("110,\t500,"));
    // Second record block start: 810
    assert!(output_content.contains("810,\t500,"));

    Ok(())
}

#[test]
fn test_lift_target() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("lift", "", "target_lift.psl");
    let sizes = get_path("lift", "", "chrom.sizes");
    let output = temp.path().join("target_lifted.psl");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("psl")
        .arg("lift")
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .arg("--t-sizes")
        .arg(&sizes);
    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    
    // Expected output check
    // First record: Target chr1:101-200 (+). tStart=10, tEnd=60.
    //   Lifted: chr1 (+). tStart=110, tEnd=160.
    // Second record: Target chr1:101-200 (-). tStart=10, tEnd=60.
    //   Lifted: chr1 (-). tStart=810, tEnd=860. (Size 1000)

    // Check first record
    // Target name: chr1
    // Target size: 1000
    // Target start: 110
    // Target end: 160
    assert!(output_content.contains("seq1\t100\t0\t50\tchr1\t1000\t110\t160"));
    
    // Check second record
    // Target name: chr1
    // Target size: 1000
    // Target start: 810
    // Target end: 860
    assert!(output_content.contains("seq1\t100\t0\t50\tchr1\t1000\t810\t860"));

    // Check tStarts
    // First: 110
    assert!(output_content.contains(",\t110,"));
    // Second: 810
    assert!(output_content.contains(",\t810,"));

    Ok(())
}
