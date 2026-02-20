use assert_cmd::cargo::cargo_bin_cmd;
use predicates::prelude::*;
use std::io::Write;

const ABCDE_NWK: &str = "((A,B),(C,D),E);";
const CATARRHINI: &str = "(((Homo,Pan),Gorilla),Pongo);";
const CATARRHINI_LABELED: &str = "(((Homo,Pan)Hominini,Gorilla)Homininae,Pongo)Hominidae;";

#[test]
fn command_prune_remove_single_leaf() -> anyhow::Result<()> {
    let mut cmd = cargo_bin_cmd!("pgr");
    cmd.arg("nwk")
        .arg("prune")
        .arg("stdin")
        .arg("-n")
        .arg("Homo")
        .write_stdin(CATARRHINI);
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("((Pan,Gorilla),Pongo);"));
    Ok(())
}

#[test]
fn command_prune_remove_multiple_leaves() -> anyhow::Result<()> {
    let mut cmd = cargo_bin_cmd!("pgr");
    cmd.arg("nwk")
        .arg("prune")
        .arg("stdin")
        .arg("-n")
        .arg("Homo")
        .arg("-n")
        .arg("Pan")
        .write_stdin(CATARRHINI);
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("(Gorilla,Pongo);"));
    Ok(())
}

#[test]
fn command_prune_remove_all_leaves_in_clade() -> anyhow::Result<()> {
    let mut cmd = cargo_bin_cmd!("pgr");
    cmd.arg("nwk")
        .arg("prune")
        .arg("stdin")
        .arg("-n")
        .arg("Homo")
        .arg("-n")
        .arg("Pan")
        .arg("-n")
        .arg("Gorilla")
        .write_stdin(CATARRHINI);
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("Pongo;"));
    Ok(())
}

#[test]
fn command_prune_remove_all_nodes_bug() -> anyhow::Result<()> {
    // Remove all nodes
    let mut cmd = cargo_bin_cmd!("pgr");
    cmd.arg("nwk")
        .arg("prune")
        .arg("stdin")
        .arg("-n")
        .arg("A")
        .arg("-n")
        .arg("B")
        .arg("-n")
        .arg("C")
        .arg("-n")
        .arg("D")
        .arg("-n")
        .arg("E")
        .write_stdin(ABCDE_NWK);
    cmd.assert().success(); // Just ensure it doesn't crash
    Ok(())
}

#[test]
fn command_prune_regex_match() -> anyhow::Result<()> {
    // Regex
    let mut cmd = cargo_bin_cmd!("pgr");
    cmd.arg("nwk")
        .arg("prune")
        .arg("stdin")
        .arg("--regex")
        .arg("^H")
        .write_stdin(CATARRHINI);
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("((Pan,Gorilla),Pongo);"));
    Ok(())
}

#[test]
fn command_prune_keep_single_leaf() -> anyhow::Result<()> {
    let mut cmd = cargo_bin_cmd!("pgr");
    cmd.arg("nwk")
        .arg("prune")
        .arg("stdin")
        .arg("-v")
        .arg("-n")
        .arg("Homo")
        .write_stdin(CATARRHINI);
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("Homo;"));
    Ok(())
}

#[test]
fn command_prune_multiple_trees() -> anyhow::Result<()> {
    let multi = format!("{}\n{}", CATARRHINI, ABCDE_NWK);
    let mut cmd = cargo_bin_cmd!("pgr");
    cmd.arg("nwk")
        .arg("prune")
        .arg("stdin")
        .arg("-n")
        .arg("Homo")
        .arg("-n")
        .arg("A")
        .write_stdin(multi);
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("((Pan,Gorilla),Pongo);"))
        .stdout(predicate::str::contains("(B,(C,D),E);"));
    Ok(())
}

#[test]
fn command_prune_file_remove_single() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    let file_path = dir.path().join("list.txt");
    {
        let mut f = std::fs::File::create(&file_path)?;
        writeln!(f, "Homo")?;
    }

    let mut cmd = cargo_bin_cmd!("pgr");
    cmd.arg("nwk")
        .arg("prune")
        .arg("stdin")
        .arg("-f")
        .arg(&file_path)
        .write_stdin(CATARRHINI);
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("((Pan,Gorilla),Pongo);"));
    Ok(())
}

#[test]
fn command_prune_file_remove_multiple() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    let file_path = dir.path().join("list.txt");
    {
        let mut f = std::fs::File::create(&file_path)?;
        writeln!(f, "Homo")?;
        writeln!(f, "Pan")?;
    }

    let mut cmd = cargo_bin_cmd!("pgr");
    cmd.arg("nwk")
        .arg("prune")
        .arg("stdin")
        .arg("-f")
        .arg(&file_path)
        .write_stdin(CATARRHINI);
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("(Gorilla,Pongo);"));
    Ok(())
}

#[test]
fn command_prune_file_remove_all() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    let file_path = dir.path().join("list.txt");
    {
        let mut f = std::fs::File::create(&file_path)?;
        writeln!(f, "Homo")?;
        writeln!(f, "Pan")?;
        writeln!(f, "Gorilla")?;
    }

    let mut cmd = cargo_bin_cmd!("pgr");
    cmd.arg("nwk")
        .arg("prune")
        .arg("stdin")
        .arg("-f")
        .arg(&file_path)
        .write_stdin(CATARRHINI);
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("Pongo;"));
    Ok(())
}

#[test]
fn command_prune_file_keep_single() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    let file_path = dir.path().join("list.txt");
    {
        let mut f = std::fs::File::create(&file_path)?;
        writeln!(f, "Homo")?;
    }

    let mut cmd = cargo_bin_cmd!("pgr");
    cmd.arg("nwk")
        .arg("prune")
        .arg("stdin")
        .arg("-v")
        .arg("-f")
        .arg(&file_path)
        .write_stdin(CATARRHINI);
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("Homo;"));
    Ok(())
}

#[test]
fn command_prune_keep_internal_node_by_label() -> anyhow::Result<()> {
    let mut cmd = cargo_bin_cmd!("pgr");
    cmd.arg("nwk")
        .arg("prune")
        .arg("stdin")
        .arg("-v")
        .arg("-n")
        .arg("Hominini")
        .write_stdin(CATARRHINI_LABELED);
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("(Homo,Pan)Hominini;"));
    Ok(())
}

#[test]
fn command_prune_keep_internal_node_by_name() -> anyhow::Result<()> {
    // Keep internal node by name, check descendants
    let mut cmd = cargo_bin_cmd!("pgr");
    cmd.arg("nwk")
        .arg("prune")
        .arg("stdin")
        .arg("-v")
        .arg("-n")
        .arg("Hominidae")
        .write_stdin(CATARRHINI_LABELED);
    // Keep Hominidae. Should keep everything under it?
    // The whole tree is Hominidae.
    cmd.assert()
        .success()
        .stdout(predicate::str::contains(CATARRHINI_LABELED));
    Ok(())
}

#[test]
fn command_prune() -> anyhow::Result<()> {
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("nwk")
        .arg("prune")
        .arg("tests/newick/catarrhini.nwk")
        .arg("-n")
        .arg("Homo")
        .arg("-n")
        .arg("Pan")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(!stdout.contains("Homo:10"));
    assert!(!stdout.contains("Gorilla:16"));
    assert!(stdout.contains("Gorilla:31"));

    Ok(())
}

#[test]
fn command_prune_invert() -> anyhow::Result<()> {
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("nwk")
        .arg("prune")
        .arg("tests/newick/catarrhini.nwk")
        .arg("-v")
        .arg("-n")
        .arg("Homo")
        .arg("-n")
        .arg("Pan")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("Homo"));
    assert!(stdout.contains("Pan"));
    assert!(!stdout.contains("Gorilla"));
    assert!(!stdout.contains("Pongo"));

    Ok(())
}
