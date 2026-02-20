use assert_cmd::cargo::cargo_bin_cmd;

#[test]
fn command_rename_basic() -> anyhow::Result<()> {
    let mut cmd = cargo_bin_cmd!("pgr");
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
    let mut cmd = cargo_bin_cmd!("pgr");
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

    let mut cmd = cargo_bin_cmd!("pgr");
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
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("nwk")
        .arg("replace")
        .arg("tests/newick/abc.nwk")
        .arg("tests/newick/abc.replace.tsv")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 1);
    assert!(stdout.contains("((Homo,Pan),Gorilla);"));

    let mut cmd = cargo_bin_cmd!("pgr");
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

    let mut cmd = cargo_bin_cmd!("pgr");
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
    let mut cmd = cargo_bin_cmd!("pgr");
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
    let mut cmd = cargo_bin_cmd!("pgr");
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
    let mut cmd = cargo_bin_cmd!("pgr");
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
    let mut cmd = cargo_bin_cmd!("pgr");
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
    let mut cmd = cargo_bin_cmd!("pgr");
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
fn command_reroot() -> anyhow::Result<()> {
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("nwk")
        .arg("reroot")
        .arg("tests/newick/catarrhini_wrong.nwk")
        .arg("-n")
        .arg("Cebus")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("(Cebus,(((Cercopithecus,(Macaca,Papio)),Simias),(Hylobates,(Pongo,(Gorilla,(Pan,Homo))))));"));

    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("nwk")
        .arg("reroot")
        .arg("tests/newick/abcde.nwk")
        .arg("-n")
        .arg("B")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("(B:0.5,(A:1,C:2)D:0.5);"));

    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("nwk")
        .arg("reroot")
        .arg("tests/newick/bs.nwk")
        .arg("-n")
        .arg("C")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("(C,(B,(A,(((D,E)86,F)93,(G,(H,I))100)100)41)61);"));

    Ok(())
}

#[test]
fn command_reroot_support() -> anyhow::Result<()> {
    // bs.nwk: (((((A,(B,C)61)41,((D,E)86,F)93)100,G)100,H,I);
    // Reroot at C with -s.
    // Labels should shift.
    // 61 (on node (B,C)) should move to node connecting (B,C) to rest?
    // Path: Root -> ... -> 41 -> 61 -> C
    // 61 moves to 41. 41 moves to ...
    // Result should show labels in new positions.
    // Just verify execution for now and check consistency.
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("nwk")
        .arg("reroot")
        .arg("tests/newick/bs.nwk")
        .arg("-n")
        .arg("C")
        .arg("-s")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // We expect the tree structure to be rerooted at C.
    // And labels shifted.
    // (C, (B, (A, ...)))
    assert!(stdout.contains("(C,"));
    assert!(stdout.contains("(B,(A,(((D,E)86,F)93,(G,(H,I)100)100)41)61)"));
    Ok(())
}

#[test]
fn command_reroot_ingroup() -> anyhow::Result<()> {
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("nwk")
        .arg("reroot")
        .arg("-l")
        .arg("tests/newick/tetrapoda.nwk")
        .arg("-n")
        .arg("Bombina")
        .arg("-n")
        .arg("Tetrao")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    println!("Output: {}", stdout);
    // Expected: ((((Procavia:0.019702...
    // Just check start and end to match test_nw_reroot_ingroup.exp content
    // (((((Procavia:0.019702,(Vulpes:0.008083,Orcinus:0.008289)84:0.008124)42:0.003924,Bradypus:0.020167)16:0.000000,(((Mesocricetus:0.011181,Tamias:0.049599)88:0.023597,Sorex:0.017660)32:0.000744,((Homo:0.004051,(Papio:0.000000,Hylobates:0.004076)42:0.000000)99:0.012677,Lepus:0.030777)67:0.007717)26:0.006246)78:0.021250,Didelphis:0.007148)71:0.0065625,(Bombina:0.269848,Tetrao:0.021544)30:0.0065625);
    assert!(stdout.contains("(((((Procavia:0.019702"));
    assert!(stdout.contains("(Bombina:0.269848,Tetrao:0.021544)30:0.0065625);"));
    Ok(())
}

#[test]
fn command_reroot_midlen() -> anyhow::Result<()> {
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("nwk")
        .arg("reroot")
        .arg("tests/newick/catarrhini.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Expected: (Pongo:15,((Gorilla:16...
    assert!(stdout.contains("(Pongo:15,((Gorilla:16"));
    Ok(())
}

#[test]
fn command_reroot_nolbl_ingrp() -> anyhow::Result<()> {
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("nwk")
        .arg("reroot")
        .arg("-l")
        .arg("tests/newick/nolbl_ingrp.nwk")
        .arg("-n")
        .arg("a")
        .arg("-n")
        .arg("b")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Input: ((a,b),c);
    // Reroot on a,b (ingroup). LCA is parent of (a,b).
    // Should result in ((a,b),c); if already rooted there?
    // newick_utils nolbl_ingrp.nw content: ((a,b),c);
    // test_nw_reroot_nolbl_ingrp.exp: ((a,b),c);
    // So it should just return the same tree or equivalent.
    assert!(stdout.contains("((a,b),c);"));
    Ok(())
}

#[test]
fn command_reroot_deroot() -> anyhow::Result<()> {
    // abcde.nwk: ((A:1,B:1)D:1,C:1)E;
    // Root is E (bifurcating, children D and C).
    // D has 2 descendants (A,B). C has 0 (leaf).
    // Deroot should splice out D.
    // Result: (A:2,B:2,C:1)E;
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("nwk")
        .arg("reroot")
        .arg("tests/newick/abcde.nwk")
        .arg("-d")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("A:2"));
    assert!(stdout.contains("B:2"));
    assert!(stdout.contains("C:1"));
    assert!(stdout.contains(")E;"));
    Ok(())
}

#[test]
fn command_reroot_lax() -> anyhow::Result<()> {
    // abcde.nwk: ((A:1,B:1)D:1,C:1)E;
    // Root is E (LCA of A and C).
    // Try to reroot on A and C. LCA is E (Root).
    // With -l, should try complement (B).
    // Reroot on B.
    // Result: (B:0.5, ...);
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("nwk")
        .arg("reroot")
        .arg("tests/newick/abcde.nwk")
        .arg("-n")
        .arg("A")
        .arg("-n")
        .arg("C")
        .arg("-l")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("B:0.5"));
    Ok(())
}

#[test]
fn command_reroot_default() -> anyhow::Result<()> {
    // abcde.nwk: ((A:1,B:1)D:1,C:1)E;
    // Max distance: A-C = 3, B-C = 3.
    // Midpoint: 1.5 from leaves.
    // From A: A->D (1) -> (0.5 on D-E edge).
    // New root is on D-E edge, 0.5 from D.
    // D branch becomes 0.5.
    // C branch becomes 0.5 (from new root to E) + 1 (E to C) = 1.5.
    // Result: ((A:1,B:1)D:0.5,C:1.5);

    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("nwk")
        .arg("reroot")
        .arg("tests/newick/abcde.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    if !stdout.contains("C:1.5") {
        println!("Output: {}", stdout);
    }
    assert!(stdout.contains("C:1.5"));
    Ok(())
}
