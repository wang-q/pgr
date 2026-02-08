use assert_cmd::Command;

#[test]
fn command_rename_basic() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("rename")
        .arg("tests/newick/catarrhini.nwk")
        .arg("-n")
        .arg("Homo")
        .arg("-r")
        .arg("Human")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("Human"));
    assert!(!stdout.contains("Homo"));
    assert!(stdout.contains("Pan")); // Others preserved

    Ok(())
}

#[test]
fn command_rename_lca() -> anyhow::Result<()> {
    // In catarrhini.nwk, Homo and Pan are children of Hominini.
    // Rename Hominini (LCA of Homo,Pan) to CladeX
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("rename")
        .arg("tests/newick/catarrhini.nwk")
        .arg("--lca")
        .arg("Homo,Pan")
        .arg("-r")
        .arg("CladeX")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("CladeX"));
    assert!(!stdout.contains("Hominini"));

    Ok(())
}

#[test]
fn command_rename_mixed() -> anyhow::Result<()> {
    // ((A,B),C);
    // Rename A -> A1.
    // Rename LCA(A,B) -> AB.

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("rename")
        .arg("stdin")
        .arg("-n")
        .arg("A")
        .arg("-r")
        .arg("A1")
        .arg("-l")
        .arg("A,B")
        .arg("-r")
        .arg("AB")
        .write_stdin("((A,B),C);")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("((A1,B)AB,C);"));

    Ok(())
}

#[test]
fn command_replace() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("replace")
        .arg("tests/newick/abc.nwk")
        .arg("tests/newick/abc.replace.tsv")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 1);
    assert!(stdout.contains("((Homo,Pan),Gorilla);"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("replace")
        .arg("tests/newick/abc.nwk")
        .arg("tests/newick/abc.replace.tsv")
        .arg("--mode")
        .arg("species")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("((A[&&NHX:S=Homo],B[&&NHX:S=Pan]),C[&&NHX:S=Gorilla]);"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("replace")
        .arg("tests/newick/abc.nwk")
        .arg("tests/newick/abc3.replace.tsv")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(
        stdout.contains("((Homo[&&NHX:color=red],Pan[&&NHX:color=red]),Gorilla[&&NHX:color=red]);")
    );

    Ok(())
}

#[test]
fn command_replace_comments() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("replace")
        .arg("tests/newick/abc.nwk")
        .arg("tests/newick/mixed_comments.replace.tsv")
        .arg("--mode")
        .arg("species")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // A -> Homo [UnstructuredComment]
    // B -> Pan [&&NHX:key=value]
    // C -> Gorilla [tag]
    assert!(stdout.contains("A[&&NHX:S=Homo:UnstructuredComment]"));
    assert!(stdout.contains("B[&&NHX:S=Pan:key=value]"));
    assert!(stdout.contains("C[&&NHX:S=Gorilla:tag]"));

    Ok(())
}

#[test]
fn command_replace_remove() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("replace")
        .arg("tests/newick/abc.nwk")
        .arg("tests/newick/replace_remove.tsv")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // A -> "" (removed)
    // B -> Pan
    // C -> C (untouched)
    // Original: ((A,B),C);
    // Expected: ((,Pan),C);
    // Note: Newick parser/writer should handle empty names as empty strings,
    // resulting in format like `(:0.1,Pan:0.1)...` or just `(,Pan)...` depending on branch lengths.
    // abc.nwk has no branch lengths.
    assert!(stdout.contains("((,Pan),C);"));

    Ok(())
}

#[test]
fn command_replace_filter() -> anyhow::Result<()> {
    // abc.nwk: ((A,B),C);
    // All are leaves.

    // 1. Skip leaves (should change nothing if all matches are leaves)
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("replace")
        .arg("tests/newick/abc.nwk")
        .arg("tests/newick/abc.replace.tsv")
        .arg("--Leaf") // Skip leaves
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // A, B, C are leaves, so they should be skipped. Output remains original (or similar).
    assert!(stdout.contains("((A,B),C);"));

    // 2. Skip internal (should change leaves)
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("replace")
        .arg("tests/newick/abc.nwk")
        .arg("tests/newick/abc.replace.tsv")
        .arg("--Internal") // Skip internal
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // A, B, C are leaves, so they should be replaced.
    assert!(stdout.contains("((Homo,Pan),Gorilla);"));

    Ok(())
}

#[test]
fn command_replace_multi() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("replace")
        .arg("tests/newick/forest.nwk")
        .arg("tests/newick/forest.replace.tsv")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Check line count (should be 5 trees)
    assert_eq!(stdout.lines().count(), 5);

    // Check replacements in the last tree (line 5)
    // Original: (Homo,(Pan,...
    // Expected: (Human,(Chimp,...
    assert!(stdout.contains("(Human,(Chimp,"));

    Ok(())
}

#[test]
fn command_subtree() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("subtree")
        .arg("tests/newick/hg38.7way.nwk")
        .arg("-n")
        .arg("Human")
        .arg("-n")
        .arg("Rhesus")
        .arg("-M")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 0);

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("subtree")
        .arg("tests/newick/hg38.7way.nwk")
        .arg("-n")
        .arg("Human")
        .arg("-n")
        .arg("Rhesus")
        .arg("-r")
        .arg("^ch")
        .arg("-M")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 1);
    assert!(stdout.contains("((Human:0.007,Chimp:0.00684):0.027,Rhesus:0.037601):0.11;"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("subtree")
        .arg("tests/newick/hg38.7way.nwk")
        .arg("-n")
        .arg("Human")
        .arg("-n")
        .arg("Rhesus")
        .arg("-r")
        .arg("^ch")
        .arg("-M")
        .arg("-c")
        .arg("Primates")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // pgr outputs NHX style comments
    assert!(stdout.contains("Primates:0.11[&&NHX:member=3:tri=white]"));

    Ok(())
}
