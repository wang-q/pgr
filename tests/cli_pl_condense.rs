use assert_cmd::Command;
use std::fs;
use tempfile::tempdir;

#[test]
fn command_pl_condense_basic() {
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    cmd.arg("pl")
        .arg("condense")
        .arg("-t")
        .arg("tests/pipeline/strains.taxon.tsv")
        .arg("tests/pipeline/minhash.reroot.newick");
    cmd.assert().success();
}

#[test]
fn command_pl_condense_with_rank() {
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    cmd.arg("pl")
        .arg("condense")
        .arg("-t")
        .arg("tests/pipeline/strains.taxon.tsv")
        .arg("-r")
        .arg("3")
        .arg("tests/pipeline/minhash.reroot.newick");
    cmd.assert().success();
}

#[test]
fn command_pl_condense_with_output() {
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("condensed.nwk");

    let mut cmd = Command::cargo_bin("pgr").unwrap();
    cmd.arg("pl")
        .arg("condense")
        .arg("-t")
        .arg("tests/pipeline/strains.taxon.tsv")
        .arg("-o")
        .arg(&output_path)
        .arg("tests/pipeline/minhash.reroot.newick");
    cmd.assert().success();

    // Check output file exists and is not empty
    assert!(output_path.exists());
    let content = fs::read_to_string(&output_path).unwrap();
    assert!(!content.is_empty());
}

#[test]
fn command_pl_condense_with_map() {
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("condensed.nwk");

    let mut cmd = Command::cargo_bin("pgr").unwrap();
    cmd.arg("pl")
        .arg("condense")
        .arg("-t")
        .arg("tests/pipeline/strains.taxon.tsv")
        .arg("--map")
        .arg("-o")
        .arg(&output_path)
        .arg("tests/pipeline/minhash.reroot.newick");
    cmd.assert().success();

    // Check output file exists
    assert!(output_path.exists());

    // Check condensed.tsv exists in current directory
    assert!(dir.path().join("condensed.tsv").exists() || std::env::current_dir().unwrap().join("condensed.tsv").exists());

    // Clean up condensed.tsv if created in current dir
    let _ = fs::remove_file("condensed.tsv");
}

#[test]
fn command_pl_condense_genus_level() {
    // Test condensing at genus level (column 3)
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    cmd.arg("pl")
        .arg("condense")
        .arg("-t")
        .arg("tests/pipeline/strains.taxon.tsv")
        .arg("-r")
        .arg("3")
        .arg("tests/pipeline/minhash.reroot.newick");
    cmd.assert().success();
}

#[test]
fn command_pl_condense_family_level() {
    // Test condensing at family level (column 4)
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    cmd.arg("pl")
        .arg("condense")
        .arg("-t")
        .arg("tests/pipeline/strains.taxon.tsv")
        .arg("-r")
        .arg("4")
        .arg("tests/pipeline/minhash.reroot.newick");
    cmd.assert().success();
}

#[test]
fn command_pl_condense_help() {
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    cmd.arg("pl").arg("condense").arg("--help");
    cmd.assert().success();
}
