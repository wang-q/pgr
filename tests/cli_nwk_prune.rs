#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::io::Write;

const ABCDE_NWK: &str = "((A,B),(C,D),E);";
const CATARRHINI: &str = "(((Homo,Pan),Gorilla),Pongo);";
const CATARRHINI_LABELED: &str = "(((Homo,Pan)Hominini,Gorilla)Homininae,Pongo)Hominidae;";

#[test]
fn command_prune_remove_single_leaf() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "prune", "stdin", "-n", "Homo"])
        .stdin(CATARRHINI)
        .run();

    assert!(stdout.contains("((Pan,Gorilla),Pongo);"));
}

#[test]
fn command_prune_remove_multiple_leaves() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "prune", "stdin", "-n", "Homo", "-n", "Pan"])
        .stdin(CATARRHINI)
        .run();

    assert!(stdout.contains("(Gorilla,Pongo);"));
}

#[test]
fn command_prune_remove_all_leaves_in_clade() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk", "prune", "stdin", "-n", "Homo", "-n", "Pan", "-n", "Gorilla",
        ])
        .stdin(CATARRHINI)
        .run();

    assert!(stdout.contains("Pongo;"));
}

#[test]
fn command_prune_remove_all_nodes_bug() {
    // Remove all nodes
    PgrCmd::new()
        .args(&[
            "nwk", "prune", "stdin", "-n", "A", "-n", "B", "-n", "C", "-n", "D", "-n", "E",
        ])
        .stdin(ABCDE_NWK)
        .run(); // Just ensure it doesn't crash
}

#[test]
fn command_prune_regex_match() {
    // Regex
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "prune", "stdin", "--regex", "^H"])
        .stdin(CATARRHINI)
        .run();

    assert!(stdout.contains("((Pan,Gorilla),Pongo);"));
}

#[test]
fn command_prune_keep_single_leaf() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "prune", "stdin", "-n", "Homo", "--invert"])
        .stdin(CATARRHINI)
        .run();

    assert!(stdout.contains("Homo;"));
}

#[test]
fn command_prune_multiple_trees() {
    let multi = format!("{}\n{}", CATARRHINI, ABCDE_NWK);
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "prune", "stdin", "-n", "Homo", "-n", "A"])
        .stdin(multi)
        .run();

    assert!(stdout.contains("((Pan,Gorilla),Pongo);"));
    assert!(stdout.contains("(B,(C,D),E);"));
}

#[test]
fn command_prune_file_remove_single() {
    let mut file = tempfile::Builder::new().suffix(".txt").tempfile().unwrap();
    writeln!(file, "Homo").unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "prune", "stdin", "-f", file.path().to_str().unwrap()])
        .stdin(CATARRHINI)
        .run();

    assert!(stdout.contains("((Pan,Gorilla),Pongo);"));
}

#[test]
fn command_prune_file_remove_multiple() {
    let mut file = tempfile::Builder::new().suffix(".txt").tempfile().unwrap();
    writeln!(file, "Homo").unwrap();
    writeln!(file, "Pan").unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "prune", "stdin", "-f", file.path().to_str().unwrap()])
        .stdin(CATARRHINI)
        .run();

    assert!(stdout.contains("(Gorilla,Pongo);"));
}

#[test]
fn command_prune_file_remove_all() {
    let mut file = tempfile::Builder::new().suffix(".txt").tempfile().unwrap();
    writeln!(file, "Homo").unwrap();
    writeln!(file, "Pan").unwrap();
    writeln!(file, "Gorilla").unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "prune", "stdin", "-f", file.path().to_str().unwrap()])
        .stdin(CATARRHINI)
        .run();

    assert!(stdout.contains("Pongo;"));
}

#[test]
fn command_prune_file_keep_single() {
    let mut file = tempfile::Builder::new().suffix(".txt").tempfile().unwrap();
    writeln!(file, "Homo").unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "prune",
            "stdin",
            "-f",
            file.path().to_str().unwrap(),
            "--invert",
        ])
        .stdin(CATARRHINI)
        .run();

    assert!(stdout.contains("Homo;"));
}

#[test]
fn command_prune_keep_internal_node_by_label() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "prune", "stdin", "--invert", "-n", "Hominini"])
        .stdin(CATARRHINI_LABELED)
        .run();

    assert!(stdout.contains("(Homo,Pan)Hominini;"));
}

#[test]
fn command_prune_keep_internal_node_by_name() {
    // Keep internal node by name, check descendants
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "prune", "stdin", "--invert", "-n", "Hominidae"])
        .stdin(CATARRHINI_LABELED)
        .run();

    // Keep Hominidae. Should keep everything under it?
    // The whole tree is Hominidae.
    assert!(stdout.contains(CATARRHINI_LABELED));
}

#[test]
fn command_prune() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "prune",
            "tests/newick/catarrhini.nwk",
            "-n",
            "Homo",
            "-n",
            "Pan",
        ])
        .run();

    assert!(!stdout.contains("Homo:10"));
    assert!(!stdout.contains("Gorilla:16"));
    assert!(stdout.contains("Gorilla:31"));
}

#[test]
fn command_prune_invert() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "prune",
            "tests/newick/catarrhini.nwk",
            "--invert",
            "-n",
            "Homo",
            "-n",
            "Pan",
        ])
        .run();

    assert!(stdout.contains("Homo"));
    assert!(stdout.contains("Pan"));
    assert!(!stdout.contains("Gorilla"));
    assert!(!stdout.contains("Pongo"));
}
