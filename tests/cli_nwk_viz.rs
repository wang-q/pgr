use assert_cmd::Command;

#[test]
fn command_indent() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("indent")
        .arg("tests/newick/hg38.7way.nwk")
        .arg("--text")
        .arg(".   ")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 19);
    assert!(stdout.contains(".   .   Human:"));
    assert!(stdout.contains("\n.   Opossum:"));

    Ok(())
}

#[test]
fn command_indent_compact() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("indent")
        .arg("tests/newick/catarrhini.nwk")
        .arg("--compact")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 1);
    assert_eq!(stdout.trim().lines().count(), 1); // Ensure only one line after trim
    assert!(stdout.contains("Gorilla"));

    Ok(())
}

#[test]
fn command_indent_simple() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("indent")
        .arg("tests/newick/catarrhini_wrong.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("  Homo,"));
    assert!(stdout.contains("      Gorilla,"));
    assert_eq!(stdout.lines().count(), 28);

    Ok(())
}

#[test]
fn command_indent_optt() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("indent")
        .arg("tests/newick/catarrhini_wrong.nwk")
        .arg("--text")
        .arg(".  ")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains(".  Homo,"));
    assert!(stdout.contains(".  .  .  Gorilla,"));
    assert_eq!(stdout.lines().count(), 28);

    Ok(())
}

#[test]
fn command_indent_multiple_optc() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("indent")
        .arg("tests/newick/forest_ind.nwk")
        .arg("--compact")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    let lines: Vec<&str> = stdout.trim().lines().collect();
    assert_eq!(lines.len(), 5);
    assert!(lines[0].starts_with("(Pandion,"));
    assert!(lines[4].starts_with("(Homo,"));

    Ok(())
}

#[test]
fn command_indent_stdin() -> anyhow::Result<()> {
    // 1. Default indentation
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("indent")
        .arg("stdin")
        .write_stdin("((A,B),C);")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Should have newlines and spaces (default 2 spaces)
    assert!(stdout.contains("  A"));
    assert!(stdout.contains("  B"));
    assert!(stdout.contains("C"));

    Ok(())
}
