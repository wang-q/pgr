#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

#[test]
fn command_indent() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "indent",
            "tests/newick/hg38.7way.nwk",
            "--text",
            ".   ",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 19);
    assert!(stdout.contains(".   .   Human:"));
    assert!(stdout.contains("\n.   Opossum:"));
}

#[test]
fn command_indent_compact() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "indent", "tests/newick/catarrhini.nwk", "--compact"])
        .run();

    assert_eq!(stdout.lines().count(), 1);
    assert_eq!(stdout.trim().lines().count(), 1); // Ensure only one line after trim
    assert!(stdout.contains("Gorilla"));
}

#[test]
fn command_indent_simple() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "indent", "tests/newick/catarrhini_wrong.nwk"])
        .run();

    assert!(stdout.contains("  Homo,"));
    assert!(stdout.contains("      Gorilla,"));
    assert_eq!(stdout.lines().count(), 28);
}

#[test]
fn command_indent_optt() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "indent",
            "tests/newick/catarrhini_wrong.nwk",
            "--text",
            ".  ",
        ])
        .run();

    assert!(stdout.contains(".  Homo,"));
    assert!(stdout.contains(".  .  .  Gorilla,"));
    assert_eq!(stdout.lines().count(), 28);
}

#[test]
fn command_indent_multiple_optc() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "indent", "tests/newick/forest_ind.nwk", "--compact"])
        .run();

    let lines: Vec<&str> = stdout.trim().lines().collect();
    assert_eq!(lines.len(), 5);
    assert!(lines[0].starts_with("(Pandion,"));
    assert!(lines[4].starts_with("(Homo,"));
}

#[test]
fn command_indent_stdin() {
    // 1. Default indentation
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "indent", "stdin"])
        .stdin("((A,B),C);")
        .run();

    // Should have newlines and spaces (default 2 spaces)
    assert!(stdout.contains("  A"));
    assert!(stdout.contains("  B"));
    assert!(stdout.contains("C"));
}

#[test]
fn command_indent_special_chars() {
    // 1. Plus/Minus in labels (plusminus.nw)
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "indent", "tests/newick/plusminus.nwk"])
        .run();

    // pgr should output it, likely quoted if it contains special chars that require quoting.
    // + is not strictly a special char in Newick (unlike (),:;), but some parsers might quote it.
    // pgr quote_label: "(),:;[] \t\n".contains(c) -> quotes.
    // + is NOT in that list. So it should be unquoted.
    assert!(stdout.contains("HRV-A+A2"));

    // 2. Slash and Space (slash_and_space.nw)
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "indent", "tests/newick/slash_and_space.nwk"])
        .run();

    // Label: B/Washington/05/2009 gi_255529494 gb_GQ451489
    // Contains space, so pgr WILL quote it.
    // newick_utils might not quote it if it's lax, but pgr is safer.
    // We just check if the text is present.
    assert!(stdout.contains("B/Washington/05/2009 gi_255529494 gb_GQ451489"));
    // Check if it is quoted
    assert!(stdout.contains("'B/Washington/05/2009 gi_255529494 gb_GQ451489'"));
}

#[test]
fn command_indent_multiple_trees() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "indent", "tests/newick/forest.nwk"])
        .run();

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
}

#[test]
fn command_comment() {
    // This test involves piping output from one pgr command to another.
    // 1. Run pgr nwk comment ... --color green
    let (color_stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "comment",
            "tests/newick/abc.nwk",
            "-n",
            "A",
            "-n",
            "C",
            "--color",
            "green",
        ])
        .run();

    // 2. Run pgr nwk comment stdin ... --dot with input from step 1
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "comment", "stdin", "-l", "A,B", "--dot"])
        .stdin(color_stdout)
        .run();

    assert_eq!(
        stdout.trim(),
        "((A[&&NHX:color=green],B)[&&NHX:dot=black],C[&&NHX:color=green]);"
    );
}

#[test]
fn command_comment_remove() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "comment",
            "tests/newick/abc.comment.nwk",
            "--remove",
            "color=",
        ])
        .run();

    assert_eq!(stdout.trim(), "((A,B)[&&NHX:dot=black],C);");
}

#[test]
fn command_to_dot() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "to-dot", "tests/newick/catarrhini.nwk"])
        .run();

    assert!(stdout.contains("digraph Tree {"));
    assert!(stdout.contains("node [shape=box];"));
    assert!(stdout.contains("Hominidae"));
    assert!(stdout.contains("->"));
}

#[test]
fn command_to_forest() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "to-forest", "tests/newick/catarrhini.nwk"])
        .run();

    assert!(stdout.contains("[,, tier="));
    assert!(stdout.contains("Hominidae"));
    assert!(stdout.contains("{Homo}"));
}

#[test]
fn command_to_forest_bl() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "to-forest", "tests/newick/catarrhini.nwk", "--bl"])
        .run();

    assert!(stdout.contains("l=")); // Should have lengths
    assert!(stdout.contains("Hominidae"));
    assert!(stdout.contains("{Homo}"));
}

#[test]
fn command_tex() {
    // 1. Default (Cladogram)
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "to-tex", "tests/newick/hg38.7way.nwk"])
        .run();

    assert!(stdout.contains(r"\documentclass"));
    assert!(stdout.contains(r"\begin{forest}"));
    assert!(stdout.contains("tier=4"));

    // 2. Phylogram (--bl)
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "to-tex", "tests/newick/hg38.7way.nwk", "--bl"])
        .run();

    assert!(stdout.contains(r"\documentclass"));
    assert!(stdout.contains("l=40mm"));
    assert!(stdout.contains("l=53mm"));
}
