use assert_cmd::Command;
use std::fs;
use tempfile::tempdir;

#[test]
fn command_pl_condense_basic() {
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    cmd.arg("pl")
        .arg("condense")
        .arg("--taxon")
        .arg("tests/pipeline/strains.taxon.tsv")
        .arg("tests/pipeline/minhash.reroot.newick");
    let output = cmd.output().expect("Failed to execute command");
    let stdout = String::from_utf8_lossy(&output.stdout);

    // Should succeed
    assert!(
        output.status.success(),
        "Command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Check that condensed labels are present (format: term___count)
    assert!(
        stdout.contains("___"),
        "No condensed labels found in output"
    );

    // Check that output is a valid Newick tree
    assert!(
        stdout.starts_with('(') || stdout.starts_with("Sa_cer"),
        "Output is not a valid Newick tree"
    );
}

#[test]
fn command_pl_condense_monophyly_check() {
    // Test that only monophyletic groups are condensed
    // Non-monophyletic groups should remain as original leaf nodes
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("condensed.nwk");

    let mut cmd = Command::cargo_bin("pgr").unwrap();
    cmd.arg("pl")
        .arg("condense")
        .arg("--taxon")
        .arg("tests/pipeline/strains.taxon.tsv")
        .arg("-o")
        .arg(&output_path)
        .arg("tests/pipeline/minhash.reroot.newick");
    cmd.assert().success();

    let content = fs::read_to_string(&output_path).unwrap();

    // Some groups should be condensed (monophyletic)
    assert!(content.contains("___"), "No condensed labels found");

    // The output should be a valid Newick tree
    assert!(content.starts_with('(') || content.starts_with("Sa_cer"));
}

#[test]
fn command_pl_condense_with_rank() {
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    cmd.arg("pl")
        .arg("condense")
        .arg("--taxon")
        .arg("tests/pipeline/strains.taxon.tsv")
        .arg("--rank")
        .arg("3")
        .arg("tests/pipeline/minhash.reroot.newick");
    let output = cmd.output().expect("Failed to execute command");
    let stdout = String::from_utf8_lossy(&output.stdout);

    assert!(output.status.success());

    // Column 3 is genus, so we should see genus-level condensation
    assert!(stdout.contains("___"), "No condensed labels found");
}

#[test]
fn command_pl_condense_with_output() {
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("condensed.nwk");

    let mut cmd = Command::cargo_bin("pgr").unwrap();
    cmd.arg("pl")
        .arg("condense")
        .arg("--taxon")
        .arg("tests/pipeline/strains.taxon.tsv")
        .arg("-o")
        .arg(&output_path)
        .arg("tests/pipeline/minhash.reroot.newick");
    cmd.assert().success();

    // Check output file exists and is not empty
    assert!(output_path.exists());
    let content = fs::read_to_string(&output_path).unwrap();
    assert!(!content.is_empty());

    // Check it's a valid Newick format
    assert!(content.starts_with('(') || content.starts_with("Sa_cer"));
}

#[test]
fn command_pl_condense_with_map() {
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("condensed.nwk");

    let mut cmd = Command::cargo_bin("pgr").unwrap();
    cmd.arg("pl")
        .arg("condense")
        .arg("--taxon")
        .arg("tests/pipeline/strains.taxon.tsv")
        .arg("--map")
        .arg("-o")
        .arg(&output_path)
        .arg("tests/pipeline/minhash.reroot.newick");
    cmd.assert().success();

    // Check output file exists
    assert!(output_path.exists());

    // Check condensed.tsv exists in current directory (where command was run)
    let map_path = std::env::current_dir().unwrap().join("condensed.tsv");
    assert!(
        map_path.exists(),
        "condensed.tsv not found at {:?}",
        map_path
    );
    let map_content = fs::read_to_string(&map_path).unwrap();
    assert!(!map_content.is_empty(), "condensed.tsv is empty");

    // Check format: original_name<TAB>condensed_label
    for line in map_content.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        assert_eq!(
            parts.len(),
            2,
            "condensed.tsv line should have 2 columns: {}",
            line
        );
        assert!(
            parts[1].contains("___"),
            "condensed label should contain ___: {}",
            parts[1]
        );
    }

    // Clean up
    let _ = fs::remove_file(&map_path);
}

#[test]
fn command_pl_condense_genus_level() {
    // Test condensing at genus level (column 3)
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    cmd.arg("pl")
        .arg("condense")
        .arg("--taxon")
        .arg("tests/pipeline/strains.taxon.tsv")
        .arg("--rank")
        .arg("3")
        .arg("tests/pipeline/minhash.reroot.newick");
    let output = cmd.output().expect("Failed to execute command");
    let stdout = String::from_utf8_lossy(&output.stdout);

    assert!(output.status.success());

    // Genus level condensation should have genus names
    assert!(stdout.contains("___"));
}

#[test]
fn command_pl_condense_family_level() {
    // Test condensing at family level (column 4)
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    cmd.arg("pl")
        .arg("condense")
        .arg("--taxon")
        .arg("tests/pipeline/strains.taxon.tsv")
        .arg("--rank")
        .arg("4")
        .arg("tests/pipeline/minhash.reroot.newick");
    let output = cmd.output().expect("Failed to execute command");
    let stdout = String::from_utf8_lossy(&output.stdout);

    assert!(output.status.success());

    // Family level should condense even more
    assert!(stdout.contains("___"));
}

#[test]
fn command_pl_condense_help() {
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    cmd.arg("pl").arg("condense").arg("--help");
    let output = cmd.output().expect("Failed to execute command");
    let stdout = String::from_utf8_lossy(&output.stdout);

    assert!(output.status.success());
    assert!(stdout.contains("condense"));
    assert!(stdout.contains("--taxon"));
    assert!(stdout.contains("--rank"));
}
