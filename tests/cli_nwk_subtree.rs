use assert_cmd::Command;

#[test]
fn command_subtree_basic() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("subtree")
        .arg("tests/newick/hg38.7way.nwk")
        .arg("-n")
        .arg("Human")
        .arg("-n")
        .arg("Rhesus")
        .arg("-M")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 0);

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("subtree")
        .arg("tests/newick/hg38.7way.nwk")
        .arg("-n")
        .arg("Human")
        .arg("-n")
        .arg("Rhesus")
        .arg("-r")
        .arg("^ch")
        .arg("-M")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 1);
    assert!(stdout.contains("((Human:0.007,Chimp:0.00684):0.027,Rhesus:0.037601):0.11;"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("subtree")
        .arg("tests/newick/hg38.7way.nwk")
        .arg("-n")
        .arg("Human")
        .arg("-n")
        .arg("Rhesus")
        .arg("-r")
        .arg("^ch")
        .arg("-M")
        .arg("-C")
        .arg("Primates")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // pgr outputs NHX style comments
    assert!(stdout.contains("Primates:0.11[&&NHX:member=3:tri=white]"));

    Ok(())
}

#[test]
fn command_subtree_context() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("subtree")
        .arg("tests/newick/hg38.7way.nwk")
        .arg("-n")
        .arg("Human")
        .arg("-n")
        .arg("Chimp")
        .arg("-c")
        .arg("1")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Context 1 should include Rhesus
    assert!(stdout.contains("Rhesus"));

    Ok(())
}

#[test]
fn command_subtree_multiple() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("subtree")
        .arg("tests/newick/catarrhini_wrong_mult.nwk")
        .arg("-n")
        .arg("Cebus")
        .arg("-n")
        .arg("Papio")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 3);
    for line in stdout.lines() {
        assert!(line.contains("(((Cercopithecus,(Macaca,Papio)),Simias),Cebus);"));
    }

    Ok(())
}

#[test]
fn command_subtree_monophyly() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("subtree")
        .arg("tests/newick/catarrhini_wrong.nwk")
        .arg("-n")
        .arg("Simias")
        .arg("-n")
        .arg("Papio")
        .arg("-n")
        .arg("Macaca")
        .arg("-n")
        .arg("Cercopithecus")
        .arg("-M")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("((Cercopithecus,(Macaca,Papio)),Simias);"));

    Ok(())
}

#[test]
fn command_subtree_monophyly_fail() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("subtree")
        .arg("tests/newick/catarrhini_wrong.nwk")
        .arg("-n")
        .arg("Simias")
        .arg("-n")
        .arg("Papio")
        .arg("-n")
        .arg("Cercopithecus")
        .arg("-M")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.trim().is_empty());

    Ok(())
}

#[test]
fn command_subtree_regex() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("subtree")
        .arg("tests/newick/HRV.nwk")
        .arg("-r")
        .arg("^HRV.*")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Only checking the structure briefly to avoid super long string matching issues
    assert!(stdout.contains("(((((HRV85_1:0.114608,(HRV89_1:0.219212,HRV1B_1:0.123339):0.076821):0.043577,"));
    // pgr formats floats without trailing zeros
    assert!(stdout.contains("HRV39_1:0.044427):0.65675,"));
    assert!(stdout.contains("):0.317738;"));

    Ok(())
}

#[test]
fn command_subtree_default() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("nwk")
        .arg("subtree")
        .arg("tests/newick/catarrhini.nwk")
        .arg("-n")
        .arg("Hominidae")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Check if it matches def_r (with branch length)
    assert!(stdout.contains("((Gorilla:16,(Pan:10,Homo:10)Hominini:10)Homininae:15,Pongo:30)Hominidae:15;"));

    Ok(())
}
