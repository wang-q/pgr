use assert_cmd::Command;
use std::process::Stdio;

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

#[test]
fn command_indent_special_chars() -> anyhow::Result<()> {
    // 1. Plus/Minus in labels (plusminus.nw)
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("indent")
        .arg("tests/newick/plusminus.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // pgr should output it, likely quoted if it contains special chars that require quoting.
    // + is not strictly a special char in Newick (unlike (),:;), but some parsers might quote it.
    // pgr quote_label: "(),:;[] \t\n".contains(c) -> quotes.
    // + is NOT in that list. So it should be unquoted.
    assert!(stdout.contains("HRV-A+A2"));

    // 2. Slash and Space (slash_and_space.nw)
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("indent")
        .arg("tests/newick/slash_and_space.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Label: B/Washington/05/2009 gi_255529494 gb_GQ451489
    // Contains space, so pgr WILL quote it.
    // newick_utils might not quote it if it's lax, but pgr is safer.
    // We just check if the text is present.
    assert!(stdout.contains("B/Washington/05/2009 gi_255529494 gb_GQ451489"));
    // Check if it is quoted
    assert!(stdout.contains("'B/Washington/05/2009 gi_255529494 gb_GQ451489'"));

    Ok(())
}

#[test]
fn command_indent_multiple_trees() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("indent")
        .arg("tests/newick/forest.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // forest.nwk contains multiple trees (5 lines).
    // pgr should output all of them.
    // Verify specific labels from different trees to ensure all are processed.
    assert!(stdout.contains("Pandion")); // From tree 1
    assert!(stdout.contains("Diomedea")); // From tree 2
    assert!(stdout.contains("Ticodendraceae")); // From tree 3
    assert!(stdout.contains("Gorilla")); // From tree 4
    assert!(stdout.contains("Cebus")); // From tree 5

    // Verify we have at least 5 semicolons (one per tree)
    assert!(stdout.matches(';').count() >= 5);

    Ok(())
}

#[test]
fn command_comment() -> anyhow::Result<()> {
    let bin = assert_cmd::cargo::cargo_bin("pgr");
    let mut cmd_color = std::process::Command::new(&bin);
    let mut child_color = cmd_color
        .arg("nwk")
        .arg("comment")
        .arg("tests/newick/abc.nwk")
        .arg("-n")
        .arg("A")
        .arg("-n")
        .arg("C")
        .arg("--color")
        .arg("green")
        .stdout(Stdio::piped())
        .spawn()?;

    let mut cmd_dot = std::process::Command::new(&bin);
    let output = cmd_dot
        .arg("nwk")
        .arg("comment")
        .arg("stdin")
        .arg("-l")
        .arg("A,B")
        .arg("--dot")
        .stdin(Stdio::from(child_color.stdout.take().unwrap()))
        .stdout(Stdio::piped())
        .spawn()?
        .wait_with_output()?;

    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(
        stdout.lines().next().unwrap(),
        "((A[&&NHX:color=green],B)[&&NHX:dot=black],C[&&NHX:color=green]);"
    );

    Ok(())
}

#[test]
fn command_comment_remove() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("comment")
        .arg("tests/newick/abc.comment.nwk")
        .arg("--remove")
        .arg("color=")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(
        stdout.lines().next().unwrap(),
        "((A,B)[&&NHX:dot=black],C);"
    );

    Ok(())
}

#[test]
fn command_to_dot() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("to-dot")
        .arg("tests/newick/catarrhini.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("digraph Tree {"));
    assert!(stdout.contains("node [shape=box];"));
    assert!(stdout.contains("Hominidae"));
    assert!(stdout.contains("->"));

    Ok(())
}

#[test]
fn command_to_forest() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("to-forest")
        .arg("tests/newick/catarrhini.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("[,, tier="));
    assert!(stdout.contains("Hominidae"));
    assert!(stdout.contains("{Homo}"));

    Ok(())
}

#[test]
fn command_to_forest_bl() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("to-forest")
        .arg("tests/newick/catarrhini.nwk")
        .arg("--bl")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("l=")); // Should have lengths
    assert!(stdout.contains("Hominidae"));
    assert!(stdout.contains("{Homo}"));

    Ok(())
}

#[test]
fn command_tex() -> anyhow::Result<()> {
    // 1. Default (Cladogram)
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("to-tex")
        .arg("tests/newick/hg38.7way.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains(r"\documentclass"));
    assert!(stdout.contains(r"\begin{forest}"));
    assert!(stdout.contains("tier=4"));

    // 2. Phylogram (--bl)
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("to-tex")
        .arg("tests/newick/hg38.7way.nwk")
        .arg("--bl")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains(r"\documentclass"));
    assert!(stdout.contains("l=40mm"));
    assert!(stdout.contains("l=53mm"));

    Ok(())
}
