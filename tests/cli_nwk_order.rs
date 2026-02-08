use assert_cmd::Command;

#[test]
fn command_order_basic() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("tests/newick/abc.nwk")
        .arg("--nd")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("(C,(A,B));"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("tests/newick/abc.nwk")
        .arg("--ndr")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("((A,B),C);"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("tests/newick/abc.nwk")
        .arg("--an")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("((A,B),C);"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("tests/newick/abc.nwk")
        .arg("--anr")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("(C,(B,A));"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("tests/newick/abc.nwk")
        .arg("--anr")
        .arg("--ndr")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("((B,A),C);"));

    Ok(())
}

#[test]
fn command_order_list() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("tests/newick/abcde.nwk")
        .arg("--list")
        .arg("tests/newick/abcde.list")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("(C:1,(B:1,A:1)D:1)E;"));

    Ok(())
}

#[test]
fn command_order_unnamed() -> anyhow::Result<()> {
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

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("stdin")
        .arg("--an")
        .write_stdin("((C,D),(A,B));")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("((A,B),(C,D));"));

    Ok(())
}

#[test]
fn command_order_species() -> anyhow::Result<()> {
    // Create a temporary directory for testing
    let tempdir = tempfile::tempdir()?;
    let temp_path = tempdir.path();

    std::fs::copy("tests/newick/species.nwk", temp_path.join("species.nwk"))?;

    // Generate a list of labels from the tree
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("nwk")
        .arg("label")
        .arg("species.nwk")
        .arg("-o")
        .arg("species.list")
        .current_dir(temp_path)
        .output()?;

    // Order the tree using the generated list
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("species.nwk")
        .arg("--list")
        .arg("species.list")
        .current_dir(temp_path)
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Compare the ordered tree with the original one
    // They should be identical as the list was generated from the original order
    let original = std::fs::read_to_string("tests/newick/species.nwk")?;
    assert_eq!(stdout.trim(), original.trim());

    // gene tree
    std::fs::copy("tests/newick/pmxc.nwk", temp_path.join("pmxc.nwk"))?;

    // Order pmxc.nwk using the generated list
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("pmxc.nwk")
        .arg("--list")
        .arg("species.list")
        .current_dir(temp_path)
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Read the original pmxc.nwk file
    let original = std::fs::read_to_string("tests/newick/pmxc.nwk")?;

    // The ordered tree should be different from the original one
    assert_ne!(stdout.trim(), original.trim());

    Ok(())
}

#[test]
fn command_order_default_catarrhini() -> anyhow::Result<()> {
    // def:catarrhini.nw
    // Expected: test_nw_order_def.exp
    let expected = "(((Cercopithecus:10,(Macaca:10,Papio:10):20)Cercopithecinae:25,(Colobus:7,Simias:10)Colobinae:5)Cercopithecidae:10,(((Gorilla:16,(Homo:10,Pan:10)Hominini:10)Homininae:15,Pongo:30)Hominidae:15,Hylobates:20):10);";

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("tests/newick/catarrhini.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.trim(), expected);
    Ok(())
}

#[test]
fn command_order_multiple_trees() -> anyhow::Result<()> {
    // mult:catarrhini_wrong_mult.nw
    // Expected: test_nw_order_mult.exp
    let expected = r#"((((((Cebus,((Cercopithecus,(Macaca,Papio)),Simias)),Hylobates),Pongo),Gorilla),Pan),Homo);
((((((Cebus,((Cercopithecus,(Macaca,Papio)),Simias)),Hylobates),Pongo),Gorilla),Pan),Homo);
((((((Cebus,((Cercopithecus,(Macaca,Papio)),Simias)),Hylobates),Pongo),Gorilla),Pan),Homo);"#;

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("tests/newick/catarrhini_wrong_mult.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Normalize newlines
    let stdout = stdout.replace("\r\n", "\n");
    assert_eq!(stdout.trim(), expected);
    Ok(())
}

#[test]
fn command_order_descendants_tetrapoda() -> anyhow::Result<()> {
    // num: -c n tetrapoda.nw
    // Expected: test_nw_order_num.exp (adjusted for Rust float formatting)
    let expected = "(Tetrao:0.015266,(Bombina:0.269848,(Didelphis:0.007148,((Bradypus:0.020167,(Procavia:0.019702,(Vulpes:0.008083,Orcinus:0.008289)84:0.008124)42:0.003924)16:0,((Sorex:0.01766,(Mesocricetus:0.011181,Tamias:0.049599)88:0.023597)32:0.000744,(Lepus:0.030777,(Homo:0.004051,(Papio:0,Hylobates:0.004076)42:0)99:0.012677)67:0.007717)26:0.006246)78:0.02125)71:0.013125)30:0.006278)100;";

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("tests/newick/tetrapoda.nwk")
        .arg("--nd")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.trim(), expected);
    Ok(())
}

#[test]
fn command_order_deladderize_verify() -> anyhow::Result<()> {
    // dl: -c d top_heavy_ladder.nw
    // Expected: test_nw_order_dl.exp
    let expected = "(Petromyzon,((Xenopus,((Equus,Homo)Mammalia,Columba)Amniota)Tetrapoda,Carcharodon)Gnathostomata)Vertebrata;";

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("order")
        .arg("tests/newick/top_heavy_ladder.nwk")
        .arg("--deladderize")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.trim(), expected);
    Ok(())
}
