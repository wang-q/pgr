#[macro_use]
#[path = "common/mod.rs"]
mod common;
use common::PgrCmd;

#[test]
fn command_order_basic() {
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "order", "tests/newick/abc.nwk", "--nd"])
        .run();

    assert!(stdout.contains("(C,(A,B));"));

    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "order", "tests/newick/abc.nwk", "--ndr"])
        .run();

    assert!(stdout.contains("((A,B),C);"));

    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "order", "tests/newick/abc.nwk", "--an"])
        .run();

    assert!(stdout.contains("((A,B),C);"));

    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "order", "tests/newick/abc.nwk", "--anr"])
        .run();

    assert!(stdout.contains("(C,(B,A));"));

    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "order", "tests/newick/abc.nwk", "--anr", "--ndr"])
        .run();

    assert!(stdout.contains("((B,A),C);"));
}

#[test]
fn command_order_list() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "order",
            "tests/newick/abcde.nwk",
            "--list",
            "tests/newick/abcde.list",
        ])
        .run();

    assert!(stdout.contains("(C:1,(B:1,A:1)D:1)E;"));
}

#[test]
fn command_order_unnamed() {
    // Test case where internal nodes are unnamed.
    // ((C,D),(A,B));
    // Without recursive label resolution:
    // (C,D) -> (C,D)
    // (A,B) -> (A,B)
    // Root -> ((C,D),(A,B)) because "" == ""
    //
    // With recursive resolution:
    // (C,D) -> rep "C"
    // (A,B) -> rep "A"
    // Root -> compares "C" vs "A", should be ((A,B),(C,D))

    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "order", "stdin", "--an"])
        .stdin("((C,D),(A,B));")
        .run();

    // With --an, it should be sorted
    assert!(stdout.contains("((A,B),(C,D));"));
}

#[test]
fn command_order_species() {
    // Create a temporary directory for testing
    let tempdir = tempfile::tempdir().unwrap();
    let temp_path = tempdir.path();

    std::fs::copy("tests/newick/species.nwk", temp_path.join("species.nwk")).unwrap();

    // Generate a list of labels from the tree
    PgrCmd::new()
        .args(&["nwk", "label", "species.nwk", "-o", "species.list"])
        .current_dir(temp_path)
        .assert()
        .success();

    // Order the tree using the generated list
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "order", "species.nwk", "--list", "species.list"])
        .current_dir(temp_path)
        .run();

    // Compare the ordered tree with the original one
    // They should be identical as the list was generated from the original order
    let original = std::fs::read_to_string("tests/newick/species.nwk").unwrap();
    assert_eq!(stdout.trim(), original.trim());

    // gene tree
    std::fs::copy("tests/newick/pmxc.nwk", temp_path.join("pmxc.nwk")).unwrap();

    // Order pmxc.nwk using the generated list
    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "order", "pmxc.nwk", "--list", "species.list"])
        .current_dir(temp_path)
        .run();

    // Read the original pmxc.nwk file
    let original = std::fs::read_to_string("tests/newick/pmxc.nwk").unwrap();

    // The ordered tree should be different from the original one
    assert_ne!(stdout.trim(), original.trim());
}

#[test]
fn command_order_default_catarrhini() {
    // def:catarrhini.nw
    // Expected: test_nw_order_def.exp
    let expected = "(((Cercopithecus:10,(Macaca:10,Papio:10):20)Cercopithecinae:25,(Colobus:7,Simias:10)Colobinae:5)Cercopithecidae:10,(((Gorilla:16,(Homo:10,Pan:10)Hominini:10)Homininae:15,Pongo:30)Hominidae:15,Hylobates:20):10);";

    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "order", "tests/newick/catarrhini.nwk"])
        .run();

    assert_eq!(stdout.trim(), expected);
}

#[test]
fn command_order_multiple_trees() {
    // mult:catarrhini_wrong_mult.nw
    // Expected: test_nw_order_mult.exp
    let expected = r#"((((((Cebus,((Cercopithecus,(Macaca,Papio)),Simias)),Hylobates),Pongo),Gorilla),Pan),Homo);
((((((Cebus,((Cercopithecus,(Macaca,Papio)),Simias)),Hylobates),Pongo),Gorilla),Pan),Homo);
((((((Cebus,((Cercopithecus,(Macaca,Papio)),Simias)),Hylobates),Pongo),Gorilla),Pan),Homo);"#;

    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "order", "tests/newick/catarrhini_wrong_mult.nwk"])
        .run();

    // Normalize newlines
    let stdout = stdout.replace("\r\n", "\n");
    assert_eq!(stdout.trim(), expected);
}

#[test]
fn command_order_descendants_tetrapoda() {
    // num: -c n tetrapoda.nw
    // Expected: test_nw_order_num.exp (adjusted for Rust float formatting)
    let expected = "(Tetrao:0.015266,(Bombina:0.269848,(Didelphis:0.007148,((Bradypus:0.020167,(Procavia:0.019702,(Vulpes:0.008083,Orcinus:0.008289)84:0.008124)42:0.003924)16:0,((Sorex:0.01766,(Mesocricetus:0.011181,Tamias:0.049599)88:0.023597)32:0.000744,(Lepus:0.030777,(Homo:0.004051,(Papio:0,Hylobates:0.004076)42:0)99:0.012677)67:0.007717)26:0.006246)78:0.02125)71:0.013125)30:0.006278)100;";

    let (stdout, _) = PgrCmd::new()
        .args(&["nwk", "order", "tests/newick/tetrapoda.nwk", "--nd"])
        .run();

    assert_eq!(stdout.trim(), expected);
}

#[test]
fn command_order_deladderize_verify() {
    // dl: -c d top_heavy_ladder.nw
    // Expected: test_nw_order_dl.exp
    let expected = "(Petromyzon,((Xenopus,((Equus,Homo)Mammalia,Columba)Amniota)Tetrapoda,Carcharodon)Gnathostomata)Vertebrata;";

    let (stdout, _) = PgrCmd::new()
        .args(&[
            "nwk",
            "order",
            "tests/newick/top_heavy_ladder.nwk",
            "--deladderize",
        ])
        .run();

    assert_eq!(stdout.trim(), expected);
}
