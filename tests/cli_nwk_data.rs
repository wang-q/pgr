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
    assert!(stdout.contains("internal labels\t5"));

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

    assert!(stdout.contains("phylogram\t19\t10\t9\t10\t5"));

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
