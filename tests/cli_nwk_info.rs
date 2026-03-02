#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::io::Write;
use tempfile::Builder;

// ================================================================================================
// pgr nwk stat
// ================================================================================================

#[test]
fn command_stat_basic() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "stat", "tests/newick/hg38.7way.nwk"])
        .run();

    assert_eq!(stdout.lines().count(), 10);
    assert!(stdout.contains("leaf labels\t7"));
    assert!(stdout.contains("rooted\tYes"));
    assert!(stdout.contains("cherries\t"));
    assert!(stdout.contains("sackin\t"));
    assert!(stdout.contains("colless\t"));
}

#[test]
fn command_stat_catarrhini() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "stat", "tests/newick/catarrhini.nwk"])
        .run();

    assert!(stdout.contains("Type\tphylogram"));
    assert!(stdout.contains("nodes\t19"));
    assert!(stdout.contains("leaves\t10"));
    assert!(stdout.contains("rooted\tYes"));
    assert!(stdout.contains("dichotomies\t9"));
    assert!(stdout.contains("leaf labels\t10"));
    assert!(stdout.contains("internal labels\t6"));
    assert!(stdout.contains("cherries\t"));
    assert!(stdout.contains("sackin\t"));
    assert!(stdout.contains("colless\t"));
}

#[test]
fn command_stat_style_line() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "stat",
            "tests/newick/catarrhini.nwk",
            "--style",
            "line",
        ])
        .run();

    assert!(stdout.contains("phylogram\t19\t10\tYes\t9\t10\t6"));
    // Header check
    assert!(stdout.contains(
        "Type\tnodes\tleaves\trooted\tdichotomies\tleaf labels\tinternal labels\tcherries\tsackin\tcolless"
    ));
}

#[test]
fn command_stat_forest() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "stat", "tests/newick/forest.nwk", "--style", "line"])
        .run();

    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 6);

    // Header
    assert!(lines[0].contains(
        "Type\tnodes\tleaves\trooted\tdichotomies\tleaf labels\tinternal labels\tcherries\tsackin\tcolless"
    ));

    // Tree 1: Cladogram, 18 nodes, 11 leaves, No rooted, 5 dichotomies, 11 leaf labels, 0 inner labels
    // 5 cherries (visual inspection of forest.nwk or just accept changes)
    assert!(lines[1].contains("cladogram\t18\t11\tNo\t5\t11\t0"));

    // Tree 2: Cladogram, 13 nodes, 8 leaves, No rooted, 3 dichotomies, 8 leaf labels, 0 inner labels
    assert!(lines[2].contains("cladogram\t13\t8\tNo\t3\t8\t0"));

    // Tree 3: Phylogram, 10 nodes, 6 leaves, No rooted, 3 dichotomies, 6 leaf labels, 0 inner labels
    assert!(lines[3].contains("phylogram\t10\t6\tNo\t3\t6\t0"));

    // Tree 4: Phylogram, 19 nodes, 10 leaves, 9 dichotomies, 10 leaf labels, 6 inner labels
    assert!(lines[4].contains("phylogram\t19\t10\tYes\t9\t10\t6"));

    // Tree 5: Cladogram, 19 nodes, 10 leaves, 9 dichotomies, 10 leaf labels, 0 inner labels
    assert!(lines[5].contains("cladogram\t19\t10\tYes\t9\t10\t0"));
}

#[test]
fn command_stat_stdin() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "stat", "stdin"])
        .stdin("((A:1,B:1):1,C:2);")
        .run();

    assert!(stdout.contains("nodes\t5"));
    assert!(stdout.contains("leaves\t3"));
    assert!(stdout.contains("leaf labels\t3"));
}

#[test]
fn command_stat_multi_tree_stdin() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "stat", "stdin"])
        .stdin("(A,B)C;(D,E)F;")
        .run();

    // Should appear twice (once for each tree)
    assert_eq!(stdout.matches("nodes\t3").count(), 2);
    assert_eq!(stdout.matches("leaves\t2").count(), 2);
}

#[test]
fn command_stat_outfile() {
    let temp_file = Builder::new().suffix(".tsv").tempfile().unwrap();
    let outfile = temp_file.path().to_str().unwrap();

    PgrCmd::new()
        .args(&["nwk", "stat", "tests/newick/catarrhini.nwk", "-o", outfile])
        .assert()
        .success();

    let content = std::fs::read_to_string(outfile).unwrap();
    assert!(content.contains("nodes\t19"));
    assert!(content.contains("leaves\t10"));
}

// ================================================================================================
// pgr nwk label
// ================================================================================================

#[test]
fn command_label_basic() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "label", "tests/newick/hg38.7way.nwk"])
        .run();

    // hg38.7way.nwk has 7 leaves (Human, Chimp, Rhesus, Mouse, Rat, Dog, Opossum)
    // and presumably no named internal nodes.
    assert_eq!(stdout.lines().count(), 7);
    assert!(stdout.contains("Human\n"));
}

#[test]
fn command_label_leaf_only() {
    // -I: Don't print internal labels (so print leaves only)
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "label", "tests/newick/catarrhini.nwk", "-I"])
        .run();

    // catarrhini.nwk has 10 leaves and 6 internal labels
    assert_eq!(stdout.lines().count(), 10);
    assert!(stdout.contains("Homo"));
    assert!(!stdout.contains("Hominini")); // Hominini is internal
}

#[test]
fn command_label_internal_only() {
    // -L: Don't print leaf labels (so print internal only)
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "label", "tests/newick/catarrhini.nwk", "-L"])
        .run();

    assert_eq!(stdout.lines().count(), 6);
    assert!(stdout.contains("Hominini"));
    assert!(!stdout.contains("Homo"));
}

#[test]
fn command_label_empty_internal() {
    // Test on a tree with no internal labels using -L
    // hg38.7way.nwk has no internal labels
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "label", "tests/newick/hg38.7way.nwk", "-L"])
        .run();

    assert_eq!(stdout.lines().count(), 0);
}

#[test]
fn command_label_selection_node_monophyly() {
    // -n selection with -M (monophyly) and -D (descendants)
    // -n Homininae -n Pongo
    // In catarrhini.nwk, Homininae is an internal node. Pongo is a leaf (genus).
    // -D includes descendants.
    // -M checks monophyly.
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "label",
            "tests/newick/catarrhini.nwk",
            "-n",
            "Homininae",
            "-n",
            "Pongo",
            "-DM",
        ])
        .run();

    // Select Homininae and Pongo, include descendants (-D), and check monophyly (-M).
    // The output contains the 4 leaf nodes of the Hominidae clade: Gorilla, Pan, Homo, Pongo.
    assert_eq!(stdout.lines().count(), 4);
}

#[test]
fn command_label_selection_file() {
    // -f file input
    let mut temp_file = Builder::new().suffix(".txt").tempfile().unwrap();
    writeln!(temp_file, "Homo").unwrap();
    writeln!(temp_file, "Pan").unwrap();
    let list_file = temp_file.path().to_str().unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "label",
            "tests/newick/catarrhini.nwk",
            "-f",
            list_file,
        ])
        .run();

    assert_eq!(stdout.lines().count(), 2);
    assert!(stdout.contains("Homo"));
    assert!(stdout.contains("Pan"));
}

#[test]
fn command_label_regex() {
    // -r regex
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "label",
            "tests/newick/hg38.7way.nwk",
            "-r",
            "^ch", // Case insensitive by default?
        ])
        .run();

    // Should match Chimp
    assert_eq!(stdout.lines().count(), 1);
    assert!(stdout.contains("Chimp"));
}

#[test]
fn command_label_regex_case_insensitive() {
    // Verify case insensitivity explicitly
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "label",
            "tests/newick/catarrhini.nwk",
            "-r",
            "^homo", // lowercase
        ])
        .run();

    // Should match Homo
    // But NOT Hominoidea (starts with Homi, not Homo)
    assert!(stdout.contains("Homo"));
    assert!(!stdout.contains("Hominoidea"));
}

#[test]
fn command_label_columns() {
    // -c columns
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "label",
            "tests/newick/catarrhini.comment.nwk",
            "-c",
            "species",
        ])
        .run();

    // Output format: Label \t Species
    // Example: Homo \t Homo
    // We expect a tab
    assert!(stdout.contains("\tHomo\n"));
}

#[test]
fn command_label_formatting_root() {
    // --root
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "label", "tests/newick/root.nwk", "--root"])
        .run();

    assert!(stdout.trim().contains("Root"));
    assert_eq!(stdout.lines().count(), 1);
}

#[test]
fn command_label_formatting_tab() {
    // -t (tab separated on one line)
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "label", "tests/newick/catarrhini.nwk", "-t"])
        .run();

    assert_eq!(stdout.lines().count(), 1);
    assert!(stdout.contains("Homo"));
    assert!(stdout.contains('\t'));
}

#[test]
fn command_label_special_chars() {
    // Special chars (slash, space)
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "label", "tests/newick/slash_and_space.nwk"])
        .run();

    assert!(stdout.contains("B/Washington/05/2009 gi_255529494 gb_GQ451489\n"));
    assert!(stdout.contains("Swit/1562056/2009_NA\n"));
    assert!(stdout.lines().count() > 10);
}

#[test]
fn command_label_multi_tree() {
    // Multiple trees in one file, -t option
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "label", "tests/newick/forest.nwk", "-t"])
        .run();

    // forest.nwk has 5 trees, so 5 lines
    assert_eq!(stdout.lines().count(), 5);
    assert!(stdout.contains("Pandion")); // Tree 1
    assert!(stdout.contains("Diomedea")); // Tree 2
    assert!(stdout.contains("Ticodendraceae")); // Tree 3
    assert!(stdout.contains("Gorilla")); // Tree 4/5
}

// ================================================================================================
// pgr nwk distance
// ================================================================================================

#[test]
fn command_distance_root() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "distance",
            "tests/newick/catarrhini.nwk",
            "-I",
            "--mode",
            "root",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 10);
    assert!(stdout.contains("Homo\t60"));
}

#[test]
fn command_distance_parent() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "distance",
            "tests/newick/catarrhini.nwk",
            "-I",
            "--mode",
            "parent",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 10);
    assert!(stdout.contains("Homo\t10"));
}

#[test]
fn command_distance_pairwise() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "distance",
            "tests/newick/catarrhini.nwk",
            "-I",
            "--mode",
            "pairwise",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 100);
    assert!(stdout.contains("Homo\tPongo\t65"));
    assert!(stdout.contains("Pongo\tHomo\t65"));
}

#[test]
fn command_distance_lca() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "distance",
            "tests/newick/catarrhini.nwk",
            "-I",
            "--mode",
            "lca",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 100);
    assert!(stdout.contains("Homo\tPongo\t35\t30"));
    assert!(stdout.contains("Homo\tHomo\t0\t0"));
}

#[test]
fn command_distance_phylip() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "distance",
            "tests/newick/catarrhini.nwk",
            "-I",
            "--mode",
            "phylip",
        ])
        .run();

    assert!(stdout.lines().count() >= 11);
    assert!(stdout.trim().starts_with("10"));
    assert!(stdout.contains("Homo"));
    assert!(stdout.contains(" 65.000000"));
}

#[test]
fn command_distance_stdin() {
    // Topological distance (stdin input)
    let input = "((A,B)C,D)E;";
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "distance", "stdin", "--mode", "root"])
        .stdin(input)
        .run();

    assert!(stdout.contains("A\t2"));
    assert!(stdout.contains("B\t2"));
    assert!(stdout.contains("C\t1"));
    assert!(stdout.contains("D\t1"));
    assert!(stdout.contains("E\t0"));
}

#[test]
fn command_distance_reference_dist_root() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "distance", "tests/newick/dist.nwk", "--mode", "root"])
        .run();

    // Verified against newick_utils test_nw_distance_rh.exp
    assert!(stdout.contains("A\t4"));
    assert!(stdout.contains("B\t6"));
    assert!(stdout.contains("C\t3"));
    assert!(stdout.contains("D\t6"));
    assert!(stdout.contains("E\t4"));
    assert!(stdout.contains("F\t4"));
}

#[test]
fn command_distance_reference_unnamed() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "distance",
            "tests/newick/dist_meth_xpl.nwk",
            "--mode",
            "root",
        ])
        .run();

    // Verified against newick_utils test_nw_distance_nsf.exp
    // Unnamed nodes might appear as tabs with no label or internal ID.
    // Based on my experience, pgr usually outputs label if present.
    // If empty label, it outputs "\tDist".
    assert!(stdout.contains("\t3"));
    assert!(stdout.contains("\t4"));
    assert!(stdout.contains("A\t5"));
    assert!(stdout.contains("B\t4"));
}

#[test]
fn command_distance_reference_lca() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "distance", "tests/newick/dist.nwk", "--mode", "lca"])
        .run();

    // Verified against newick_utils test_nw_distance_an_2.exp (D F -> 4 2)
    // pgr output: D \t F \t 4 \t 2
    assert!(stdout.contains("D\tF\t4\t2"));
}

#[test]
fn command_distance_reference_pairwise() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "distance",
            "tests/newick/dist.nwk",
            "--mode",
            "pairwise",
        ])
        .run();

    // Verified against newick_utils test_nw_distance_pan.exp (F D E B)
    // Check F-D distance (2+4=6)
    // Check F-E distance (2+2=4)
    // Check D-B distance (6+6=12)
    assert!(stdout.contains("F\tD\t6"));
    assert!(stdout.contains("F\tE\t4"));
    assert!(stdout.contains("D\tB\t12"));
}

#[test]
fn command_distance_reference_phylip() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "distance",
            "tests/newick/dist.nwk",
            "--mode",
            "phylip",
        ])
        .run();

    // Verified against newick_utils test_nw_distance_m.exp
    // Note: pgr includes all named nodes (leaves + internal) in phylip mode.
    // dist.nwk has 6 leaves + 5 named internal nodes = 11 nodes.
    // newick_utils defaults to leaves only for matrix.
    assert!(stdout.lines().next().unwrap().trim().starts_with("11"));
    assert!(stdout.contains("A"));
    assert!(stdout.contains("B"));
    // Check for some distance values in the output
    assert!(stdout.contains("6.000000"));
    assert!(stdout.contains("7.000000"));
    assert!(stdout.contains("10.000000"));
}

#[test]
fn command_distance_float_noise() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "distance",
            "tests/newick/noise.nwk",
            "--mode",
            "pairwise",
        ])
        .run();

    // A->B should be 0.1 + 0.2 = 0.3
    // Without fix: 0.30000000000000004
    // With fix: 0.3
    assert!(stdout.contains("A\tB\t0.3\n") || stdout.contains("A\tB\t0.30\n")); // Allow formatted but clean
    assert!(!stdout.contains("0.30000000000000004"));
}

#[test]
fn command_distance_reference_parent() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "distance",
            "tests/newick/dist.nwk",
            "--mode",
            "parent",
        ])
        .run();

    // Verified against newick_utils test_nw_distance_par_all_nam.exp
    // A->g: 2, B->g: 4, g->k: 2
    // C->j: 2, D->h: 3, E->h: 1
    // h->i: 1, F->i: 2, i->j: 1
    // j->k: 1, k->None: 0
    assert!(stdout.contains("A\t2"));
    assert!(stdout.contains("B\t4"));
    assert!(stdout.contains("g\t2"));
    assert!(stdout.contains("C\t2"));
    assert!(stdout.contains("D\t3"));
    assert!(stdout.contains("E\t1"));
    assert!(stdout.contains("h\t1"));
    assert!(stdout.contains("F\t2"));
    assert!(stdout.contains("i\t1"));
    assert!(stdout.contains("j\t1"));
    assert!(stdout.contains("k\t0"));
}
