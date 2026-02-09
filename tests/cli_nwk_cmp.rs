use assert_cmd::Command;
use std::io::Write;
use tempfile::Builder;

#[test]
fn command_nwk_cmp_single_file() -> anyhow::Result<()> {
    // Create a temporary Newick file with 2 different trees
    let mut file = Builder::new().suffix(".nwk").tempfile()?;
    // Tree 1: ((A,B),(C,D)); -> Splits: {A,B} vs {C,D}
    // Tree 2: ((A,C),(B,D)); -> Splits: {A,C} vs {B,D}
    // RF distance should be 2.
    writeln!(file, "((A,B),(C,D));")?;
    writeln!(file, "((A,C),(B,D));")?;

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("nwk").arg("cmp").arg(file.path()).output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Expected output:
    // Tree1 Tree2 RF_Dist
    // 1     1     0
    // 1     2     2
    // 2     1     2
    // 2     2     0

    assert!(stdout.contains("Tree1\tTree2\tRF_Dist"));
    assert!(stdout.contains("1\t1\t0"));
    assert!(stdout.contains("1\t2\t2"));
    assert!(stdout.contains("2\t1\t2"));
    assert!(stdout.contains("2\t2\t0"));

    Ok(())
}

#[test]
fn command_nwk_cmp_two_files() -> anyhow::Result<()> {
    let mut file1 = Builder::new().suffix(".nwk").tempfile()?;
    writeln!(file1, "((A,B),(C,D));")?; // Tree 1

    let mut file2 = Builder::new().suffix(".nwk").tempfile()?;
    writeln!(file2, "((A,B),(C,D));")?; // Tree 1 (Same)
    writeln!(file2, "((A,C),(B,D));")?; // Tree 2 (Diff, RF=2)

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("cmp")
        .arg(file1.path())
        .arg(file2.path())
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Expected:
    // T1(File1) vs T1(File2) -> 0
    // T1(File1) vs T2(File2) -> 2

    assert!(stdout.contains("1\t1\t0"));
    assert!(stdout.contains("1\t2\t2"));

    Ok(())
}
