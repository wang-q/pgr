use assert_cmd::Command;
use std::fs;
use tempfile::TempDir;

#[test]
fn test_2bit_some() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test_some.fa");
    let twobit = temp.path().join("test_some.2bit");
    let list = temp.path().join("list.txt");
    let output = temp.path().join("out_some.fa");
    let output_inv = temp.path().join("out_some_inv.fa");

    // seq1: ACGT (4)
    // seq2: TGCA (4)
    // seq3: NNNN (4)
    fs::write(&input, ">seq1\nACGT\n>seq2\nTGCA\n>seq3\nNNNN\n")?;
    fs::write(&list, "seq1\nseq3\n")?;

    // Create 2bit
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("fa")
        .arg("to2bit")
        .arg(&input)
        .arg("-o")
        .arg(&twobit);
    cmd.assert().success();

    // Test some
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("some")
        .arg(&twobit)
        .arg(&list)
        .arg("-o")
        .arg(&output);
    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    assert!(output_content.contains(">seq1"));
    assert!(output_content.contains("ACGT"));
    assert!(output_content.contains(">seq3"));
    assert!(output_content.contains("NNNN"));
    assert!(!output_content.contains(">seq2"));

    // Test some invert
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("some")
        .arg(&twobit)
        .arg(&list)
        .arg("-i")
        .arg("-o")
        .arg(&output_inv);
    cmd.assert().success();

    let output_inv_content = fs::read_to_string(&output_inv)?;
    assert!(output_inv_content.contains(">seq2"));
    assert!(output_inv_content.contains("TGCA"));
    assert!(!output_inv_content.contains(">seq1"));
    assert!(!output_inv_content.contains(">seq3"));

    Ok(())
}
