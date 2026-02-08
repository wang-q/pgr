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

    assert!(stdout.contains("((Homo[&&NHX:color=red],Pan[&&NHX:color=red]),Gorilla[&&NHX:color=red]);"));

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
fn command_topo_basic() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("topo")
        .arg("tests/newick/catarrhini.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // catarrhini.nwk has lengths and comments.
    // Default topo removes lengths and comments (properties), keeps labels.
    // Original: ((((Gorilla:16,(Pan:10,Homo:10)Hominini:10)Homininae:15,Pongo:30)Hominidae:15,Hylobates:20):10,(((Macaca:10,Papio:10):20,Cercopithecus:10)Cercopithecinae:25,(Simias:10,Colobus:7)Colobinae:5)Cercopithecidae:10);
    // Expected: ((((Gorilla,(Pan,Homo)Hominini)Homininae,Pongo)Hominidae,Hylobates),(((Macaca,Papio),Cercopithecus)Cercopithecinae,(Simias,Colobus)Colobinae)Cercopithecidae);
    // Note: The root edge length is also removed.
    
    assert!(stdout.contains("((((Gorilla,(Pan,Homo)Hominini)Homininae,Pongo)Hominidae,Hylobates)"));
    assert!(!stdout.contains(":")); // No lengths

    Ok(())
}

#[test]
fn command_topo_remove_labels() -> anyhow::Result<()> {
    // Test with -I (remove internal labels) and -L (remove leaf labels)
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("topo")
        .arg("tests/newick/catarrhini.nwk")
        .arg("-I")
        .arg("-L")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Should have no labels
    assert!(stdout.contains("((((,(,))")); 
    assert!(!stdout.contains("Homo"));
    assert!(!stdout.contains("Hominini"));

    Ok(())
}

#[test]
fn command_topo_keep_bl() -> anyhow::Result<()> {
    // Test --bl (keep branch lengths)
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("topo")
        .arg("tests/newick/catarrhini.nwk")
        .arg("-I")
        .arg("-L")
        .arg("--bl")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains(":16")); // Check for specific length
    assert!(!stdout.contains("Gorilla"));

    Ok(())
}


