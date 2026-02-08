use assert_cmd::Command;
use tempfile::Builder;

#[test]
fn command_stat() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("stat")
        .arg("tests/newick/hg38.7way.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 6);
    assert!(stdout.contains("leaf labels\t7"));

    Ok(())
}

#[test]
fn command_stat_catarrhini() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("stat")
        .arg("tests/newick/catarrhini.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("Type\tphylogram"));
    assert!(stdout.contains("nodes\t19"));
    assert!(stdout.contains("leaves\t10"));
    assert!(stdout.contains("dichotomies\t9"));
    assert!(stdout.contains("leaf labels\t10"));
    assert!(stdout.contains("internal labels\t6"));

    Ok(())
}

#[test]
fn command_stat_style_line() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("stat")
        .arg("tests/newick/catarrhini.nwk")
        .arg("--style")
        .arg("line")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("phylogram\t19\t10\t9\t10\t6"));

    Ok(())
}

#[test]
fn command_stat_forest() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("stat")
        .arg("tests/newick/forest.nwk")
        .arg("--style")
        .arg("line")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 6);

    // Header
    assert!(lines[0].contains("Type\tnodes\tleaves\tdichotomies\tleaf labels\tinternal labels"));

    // Tree 1: Cladogram, 18 nodes, 11 leaves, 5 dichotomies, 11 leaf labels, 0 inner labels
    assert!(lines[1].contains("cladogram\t18\t11\t5\t11\t0"));

    // Tree 2: Cladogram, 13 nodes, 8 leaves, 3 dichotomies, 8 leaf labels, 0 inner labels
    assert!(lines[2].contains("cladogram\t13\t8\t3\t8\t0"));

    // Tree 3: Phylogram, 10 nodes, 6 leaves, 3 dichotomies, 6 leaf labels, 0 inner labels
    assert!(lines[3].contains("phylogram\t10\t6\t3\t6\t0"));

    // Tree 4: Phylogram, 19 nodes, 10 leaves, 9 dichotomies, 10 leaf labels, 6 inner labels
    assert!(lines[4].contains("phylogram\t19\t10\t9\t10\t6"));

    // Tree 5: Cladogram, 19 nodes, 10 leaves, 9 dichotomies, 10 leaf labels, 0 inner labels
    assert!(lines[5].contains("cladogram\t19\t10\t9\t10\t0"));

    Ok(())
}

#[test]
fn command_stat_multi_tree() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("stat")
        .arg("stdin")
        .write_stdin("(A,B)C;(D,E)F;")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Should appear twice (once for each tree)
    assert_eq!(stdout.matches("nodes\t3").count(), 2);
    assert_eq!(stdout.matches("leaves\t2").count(), 2);

    Ok(())
}

#[test]
fn command_stat_stdin() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("stat")
        .arg("stdin")
        .write_stdin("((A:1,B:1):1,C:2);")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("nodes\t5"));
    assert!(stdout.contains("leaves\t3"));
    assert!(stdout.contains("leaf labels\t3"));

    Ok(())
}

#[test]
fn command_stat_outfile() -> anyhow::Result<()> {
    let temp_file = Builder::new().suffix(".tsv").tempfile()?;
    let outfile = temp_file.path().to_str().unwrap();

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("nwk")
        .arg("stat")
        .arg("tests/newick/catarrhini.nwk")
        .arg("-o")
        .arg(outfile)
        .assert()
        .success();

    let content = std::fs::read_to_string(outfile)?;
    assert!(content.contains("nodes\t19"));
    assert!(content.contains("leaves\t10"));

    Ok(())
}


#[test]
fn command_label() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("label")
        .arg("tests/newick/hg38.7way.nwk")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 7);
    assert!(stdout.contains("Human\n"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("label")
        .arg("tests/newick/hg38.7way.nwk")
        .arg("-L")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 0);

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("label")
        .arg("tests/newick/hg38.7way.nwk")
        .arg("-r")
        .arg("^ch")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 1);

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("label")
        .arg("tests/newick/catarrhini.nwk")
        .arg("-n")
        .arg("Homininae")
        .arg("-n")
        .arg("Pongo")
        .arg("-DM")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 4);

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("label")
        .arg("tests/newick/catarrhini.comment.nwk")
        .arg("-c")
        .arg("species")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("\tHomo\n"));

    Ok(())
}
