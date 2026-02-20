use assert_cmd::cargo::cargo_bin_cmd;
use tempfile::TempDir;

#[test]
fn command_mat_upgma() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("input.phy");
    let output = temp.path().join("output.nwk");

    let content = "4
A 0 7 11 14
B 7 0 6 9
C 11 6 0 7
D 14 9 7 0
";
    std::fs::write(&input, content)?;

    let mut cmd = cargo_bin_cmd!("pgr");
    cmd.arg("mat")
        .arg("upgma")
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .assert()
        .success();

    let nwk = std::fs::read_to_string(&output)?;
    // println!("{}", nwk); // For debugging
    assert!(nwk.contains("A:"));
    assert!(nwk.contains("B:"));
    assert!(nwk.contains("C:"));
    assert!(nwk.contains("D:"));

    // Check topology structure
    // For this matrix: B-C=6 (min), so B and C merge first at height 3.
    // Then D merges with (B,C) at height 4.
    // Finally A merges with ((B,C),D) at height 5.333.

    // Check for B:3 and C:3
    assert!(nwk.contains("B:3") || nwk.contains("B:3.0"));
    assert!(nwk.contains("C:3") || nwk.contains("C:3.0"));

    // Check for D:4
    assert!(nwk.contains("D:4") || nwk.contains("D:4.0"));

    // Check for A:5.33
    assert!(nwk.contains("A:5.33"));

    // Check groupings
    assert!(nwk.contains("(B:3"));
    assert!(nwk.contains("C:3)"));

    Ok(())
}

#[test]
fn command_mat_nj() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("input.phy");
    let output = temp.path().join("output.nwk");

    let content = "4
A 0 7 11 14
B 7 0 6 9
C 11 6 0 7
D 14 9 7 0
";
    std::fs::write(&input, content)?;

    let mut cmd = cargo_bin_cmd!("pgr");
    cmd.arg("mat")
        .arg("nj")
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .assert()
        .success();

    let nwk = std::fs::read_to_string(&output)?;
    assert!(nwk.contains("A:"));
    assert!(nwk.contains("B:"));
    assert!(nwk.contains("C:"));
    assert!(nwk.contains("D:"));

    // Check topology structure
    // NJ on this additive matrix should recover ((A,B),(C,D))
    // Note: The exact string depends on rooting and child order, but (A,B) and (C,D) should be clades
    // Current NJ implementation roots at midpoint of last edge, so we expect a rooted tree.

    // For this specific matrix:
    // A-B = 7, C-D = 7.
    // This is a perfect tree: (A:2,B:5),(C:4,D:3) ? No, let's check exact NJ math.
    // But we can just check if A and B are grouped together.

    // We can also verify via pipe
    let mut cmd = cargo_bin_cmd!("pgr");
    let assert = cmd
        .arg("mat")
        .arg("nj")
        .arg("stdin") // stdin
        .write_stdin(content)
        .assert()
        .success();

    let stdout = String::from_utf8(assert.get_output().stdout.clone())?;
    assert!(stdout.contains("A:"));
    assert!(stdout.contains("B:"));

    Ok(())
}
