#[macro_use]
#[path = "common/mod.rs"]
mod common;
use common::PgrCmd;

#[test]
fn command_rename_basic() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "rename",
            "tests/newick/catarrhini.nwk",
            "-n",
            "Homo",
            "-r",
            "Human",
        ])
        .run();

    assert!(stdout.contains("Human"));
    assert!(!stdout.contains("Homo"));
    assert!(stdout.contains("Pan")); // Others preserved
}

#[test]
fn command_rename_lca() {
    // In catarrhini.nwk, Homo and Pan are children of Hominini.
    // Rename Hominini (LCA of Homo,Pan) to CladeX
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "rename",
            "tests/newick/catarrhini.nwk",
            "--lca",
            "Homo,Pan",
            "-r",
            "CladeX",
        ])
        .run();

    assert!(stdout.contains("CladeX"));
    assert!(!stdout.contains("Hominini"));
}

#[test]
fn command_rename_mixed() {
    // ((A,B),C);
    // Rename A -> A1.
    // Rename LCA(A,B) -> AB.

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk", "rename", "stdin", "-n", "A", "-r", "A1", "-l", "A,B", "-r", "AB",
        ])
        .stdin("((A,B),C);")
        .run();

    assert!(stdout.contains("((A1,B)AB,C);"));
}

#[test]
fn command_replace() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "replace",
            "tests/newick/abc.nwk",
            "tests/newick/abc.replace.tsv",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 1);
    assert!(stdout.contains("((Homo,Pan),Gorilla);"));

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "replace",
            "tests/newick/abc.nwk",
            "tests/newick/abc.replace.tsv",
            "--mode",
            "species",
        ])
        .run();

    assert!(stdout.contains("((A[&&NHX:S=Homo],B[&&NHX:S=Pan]),C[&&NHX:S=Gorilla]);"));

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "replace",
            "tests/newick/abc.nwk",
            "tests/newick/abc.replace.tsv",
            "--mode",
            "asis",
        ])
        .run();

    assert!(stdout.contains("((A[&&NHX:Homo],B[&&NHX:Pan]),C[&&NHX:Gorilla]);"));
}

#[test]
fn command_replace_comments() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "replace",
            "tests/newick/abc.nwk",
            "tests/newick/mixed_comments.replace.tsv",
            "--mode",
            "species",
        ])
        .run();

    // A -> Homo [UnstructuredComment]
    // B -> Pan [&&NHX:key=value]
    // C -> Gorilla [tag]
    assert!(stdout.contains("A[&&NHX:S=Homo:UnstructuredComment]"));
    assert!(stdout.contains("B[&&NHX:S=Pan:key=value]"));
    assert!(stdout.contains("C[&&NHX:S=Gorilla:tag]"));
}

#[test]
fn command_replace_remove() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "replace",
            "tests/newick/abc.nwk",
            "tests/newick/replace_remove.tsv",
        ])
        .run();

    // A -> "" (removed)
    // B -> Pan
    // C -> C (untouched)
    // Original: ((A,B),C);
    // Expected: ((,Pan),C);
    // Note: Newick parser/writer should handle empty names as empty strings,
    // resulting in format like `(:0.1,Pan:0.1)...` or just `(,Pan)...` depending on branch lengths.
    // abc.nwk has no branch lengths.
    assert!(stdout.contains("((,Pan),C);"));
}

#[test]
fn command_replace_filter() {
    // abc.nwk: ((A,B),C);
    // All are leaves.

    // 1. Skip leaves (should change nothing if all matches are leaves)
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "replace",
            "tests/newick/abc.nwk",
            "tests/newick/abc.replace.tsv",
            "--Leaf", // Skip leaves
        ])
        .run();

    // A, B, C are leaves, so they should be skipped. Output remains original (or similar).
    assert!(stdout.contains("((A,B),C);"));

    // 2. Skip internal (should change leaves)
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "replace",
            "tests/newick/abc.nwk",
            "tests/newick/abc.replace.tsv",
            "--Internal", // Skip internal
        ])
        .run();

    // A, B, C are leaves, so they should be replaced.
    assert!(stdout.contains("((Homo,Pan),Gorilla);"));
}

#[test]
fn command_replace_multi() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "replace",
            "tests/newick/forest.nwk",
            "tests/newick/forest.replace.tsv",
        ])
        .run();

    // Check line count (should be 5 trees)
    assert_eq!(stdout.lines().count(), 5);

    // Check replacements in the last tree (line 5)
    // Original: (Homo,(Pan,...
    // Expected: (Human,(Chimp,...
    assert!(stdout.contains("(Human,(Chimp,"));
}

#[test]
fn command_reroot() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "reroot",
            "tests/newick/catarrhini_wrong.nwk",
            "-n",
            "Cebus",
        ])
        .run();

    assert!(stdout.contains("(Cebus,(((Cercopithecus,(Macaca,Papio)),Simias),(Hylobates,(Pongo,(Gorilla,(Pan,Homo))))));"));

    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "reroot", "tests/newick/abcde.nwk", "-n", "B"])
        .run();

    assert!(stdout.contains("(B:0.5,(A:1,C:2)D:0.5);"));

    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "reroot", "tests/newick/bs.nwk", "-n", "C"])
        .run();

    assert!(stdout.contains("(C,(B,(A,(((D,E)86,F)93,(G,(H,I))100)100)41)61);"));
}

#[test]
fn command_reroot_support() {
    // bs.nwk: (((((A,(B,C)61)41,((D,E)86,F)93)100,G)100,H,I);
    // Reroot at C with -s.
    // Labels should shift.
    // 61 (on node (B,C)) should move to node connecting (B,C) to rest?
    // Path: Root -> ... -> 41 -> 61 -> C
    // 61 moves to 41. 41 moves to ...
    // Result should show labels in new positions.
    // Just verify execution for now and check consistency.
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "reroot", "tests/newick/bs.nwk", "-n", "C", "-s"])
        .run();

    // We expect the tree structure to be rerooted at C.
    // And labels shifted.
    // (C, (B, (A, ...)))
    assert!(stdout.contains("(C,"));
    assert!(stdout.contains("(B,(A,(((D,E)86,F)93,(G,(H,I)100)100)41)61)"));
}

#[test]
fn command_reroot_ingroup() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "reroot",
            "-l",
            "tests/newick/tetrapoda.nwk",
            "-n",
            "Bombina",
            "-n",
            "Tetrao",
        ])
        .run();

    println!("Output: {}", stdout);
    // Expected: ((((Procavia:0.019702...
    // Just check start and end to match test_nw_reroot_ingroup.exp content
    // (((((Procavia:0.019702,(Vulpes:0.008083,Orcinus:0.008289)84:0.008124)42:0.003924,Bradypus:0.020167)16:0.000000,(((Mesocricetus:0.011181,Tamias:0.049599)88:0.023597,Sorex:0.017660)32:0.000744,((Homo:0.004051,(Papio:0.000000,Hylobates:0.004076)42:0.000000)99:0.012677,Lepus:0.030777)67:0.007717)26:0.006246)78:0.021250,Didelphis:0.007148)71:0.0065625,(Bombina:0.269848,Tetrao:0.021544)30:0.0065625);
    assert!(stdout.contains("(((((Procavia:0.019702"));
    assert!(stdout.contains("(Bombina:0.269848,Tetrao:0.021544)30:0.0065625);"));
}

#[test]
fn command_reroot_midlen() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "reroot", "tests/newick/catarrhini.nwk"])
        .run();

    // Expected: (Pongo:15,((Gorilla:16...
    assert!(stdout.contains("(Pongo:15,((Gorilla:16"));
}

#[test]
fn command_reroot_nolbl_ingrp() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "reroot",
            "-l",
            "tests/newick/nolbl_ingrp.nwk",
            "-n",
            "a",
            "-n",
            "b",
        ])
        .run();

    // Input: ((a,b),c);
    // Reroot on a,b (ingroup). LCA is parent of (a,b).
    // Should result in ((a,b),c); if already rooted there?
    // newick_utils nolbl_ingrp.nw content: ((a,b),c);
    // test_nw_reroot_nolbl_ingrp.exp: ((a,b),c);
    // So it should just return the same tree or equivalent.
    assert!(stdout.contains("((a,b),c);"));
}

#[test]
fn command_reroot_deroot() {
    // abcde.nwk: ((A:1,B:1)D:1,C:1)E;
    // Root is E (bifurcating, children D and C).
    // D has 2 descendants (A,B). C has 0 (leaf).
    // Deroot should splice out D.
    // Result: (A:2,B:2,C:1)E;
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "reroot", "tests/newick/abcde.nwk", "-d"])
        .run();

    assert!(stdout.contains("(A:2,B:2,C:1)E;"));
}

#[test]
fn command_reroot_lax() -> anyhow::Result<()> {
    // abcde.nwk: ((A:1,B:1)D:1,C:1)E;
    // Root is E (LCA of A and C).
    // Try to reroot on A and C. LCA is E (Root).
    // With -l, should try complement (B).
    // Reroot on B.
    // Result: (B:0.5, ...);
    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
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

    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
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
