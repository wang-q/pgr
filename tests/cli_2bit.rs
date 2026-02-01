use assert_cmd::Command;
use std::fs;
use tempfile::TempDir;
use pgr::libs::twobit::TwoBitFile;

#[test]
fn test_2bit_to2bit() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test.fa");
    let output = temp.path().join("out.2bit");
    
    fs::write(&input, ">seq1\nACGT\n>seq2\nNNNN\n")?;
    
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("to2bit")
        .arg(&input)
        .arg("-o")
        .arg(&output);
    
    cmd.assert().success();
    
    assert!(output.exists());
    
    let mut tb = TwoBitFile::open(&output)?;
    let names = tb.get_sequence_names();
    assert_eq!(names.len(), 2);
    assert!(names.contains(&"seq1".to_string()));
    assert!(names.contains(&"seq2".to_string()));
    
    let seq1 = tb.read_sequence("seq1", None, None, false)?;
    assert_eq!(seq1, "ACGT");
    
    let seq2 = tb.read_sequence("seq2", None, None, false)?;
    assert_eq!(seq2, "NNNN");

    Ok(())
}

#[test]
fn test_2bit_to2bit_strip_version() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test_ver.fa");
    let output = temp.path().join("out_ver.2bit");
    
    fs::write(&input, ">NM_001.1\nACGT\n")?;
    
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("to2bit")
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .arg("--strip-version");
    
    cmd.assert().success();
    
    let tb = TwoBitFile::open(&output)?;
    let names = tb.get_sequence_names();
    assert_eq!(names[0], "NM_001");
    
    Ok(())
}

#[test]
fn test_2bit_to2bit_mask() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test_mask.fa");
    let output = temp.path().join("out_mask.2bit");
    
    fs::write(&input, ">seq1\nacgtACGT\n")?;
    
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("to2bit")
        .arg(&input)
        .arg("-o")
        .arg(&output);
    
    cmd.assert().success();
    
    let mut tb = TwoBitFile::open(&output)?;
    let seq_masked = tb.read_sequence("seq1", None, None, false)?;
    assert_eq!(seq_masked, "acgtACGT");
    
    let seq_unmasked = tb.read_sequence("seq1", None, None, true)?;
    assert_eq!(seq_unmasked, "ACGTACGT");

    Ok(())
}

#[test]
fn test_2bit_tofa_basic() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test.fa");
    let twobit = temp.path().join("test.2bit");
    let output = temp.path().join("out.fa");
    
    fs::write(&input, ">seq1\nACGT\n>seq2\nNNNN\n")?;
    
    // Create 2bit first
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("to2bit")
        .arg(&input)
        .arg("-o")
        .arg(&twobit);
    cmd.assert().success();
    
    // Convert back to FASTA
    let mut cmd2 = Command::cargo_bin("pgr")?;
    cmd2.arg("2bit")
        .arg("tofa")
        .arg(&twobit)
        .arg("-o")
        .arg(&output);
    cmd2.assert().success();
    
    let content = fs::read_to_string(&output)?;
    // Order might differ, but content should match.
    // >seq1
    // ACGT
    // >seq2
    // NNNN
    assert!(content.contains(">seq1\nACGT"));
    assert!(content.contains(">seq2\nNNNN"));

    Ok(())
}

#[test]
fn test_2bit_tofa_seq_range() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test.fa");
    let twobit = temp.path().join("test.2bit");
    let output = temp.path().join("out.fa");
    
    fs::write(&input, ">seq1\nACGTACGT\n")?;
    
    // Create 2bit
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("to2bit")
        .arg(&input)
        .arg("-o")
        .arg(&twobit);
    cmd.assert().success();
    
    // Extract range 1-5 (CGTA)
    let mut cmd2 = Command::cargo_bin("pgr")?;
    cmd2.arg("2bit")
        .arg("tofa")
        .arg(&twobit)
        .arg("--seq")
        .arg("seq1")
        .arg("--start")
        .arg("1")
        .arg("--end")
        .arg("5")
        .arg("-o")
        .arg(&output);
    cmd2.assert().success();
    
    let content = fs::read_to_string(&output)?;
    assert!(content.contains(">seq1:1-5\nCGTA"));

    Ok(())
}

#[test]
fn test_2bit_tofa_seq_list() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test.fa");
    let twobit = temp.path().join("test.2bit");
    let list = temp.path().join("list.txt");
    let output = temp.path().join("out.fa");
    
    fs::write(&input, ">seq1\nACGT\n>seq2\nTGCA\n")?;
    fs::write(&list, "seq2\n")?;
    
    // Create 2bit
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("to2bit")
        .arg(&input)
        .arg("-o")
        .arg(&twobit);
    cmd.assert().success();
    
    // Extract seq2
    let mut cmd2 = Command::cargo_bin("pgr")?;
    cmd2.arg("2bit")
        .arg("tofa")
        .arg(&twobit)
        .arg("--seqList")
        .arg(&list)
        .arg("-o")
        .arg(&output);
    cmd2.assert().success();
    
    let content = fs::read_to_string(&output)?;
    assert!(content.contains(">seq2\nTGCA"));
    assert!(!content.contains(">seq1"));

    Ok(())
}

#[test]
fn test_2bit_tofa_mask() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test.fa");
    let twobit = temp.path().join("test.2bit");
    let output_masked = temp.path().join("masked.fa");
    let output_unmasked = temp.path().join("unmasked.fa");
    
    fs::write(&input, ">seq1\nacgtACGT\n")?;
    
    // Create 2bit
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("to2bit")
        .arg(&input)
        .arg("-o")
        .arg(&twobit);
    cmd.assert().success();
    
    // Default (masked)
    let mut cmd2 = Command::cargo_bin("pgr")?;
    cmd2.arg("2bit")
        .arg("tofa")
        .arg(&twobit)
        .arg("-o")
        .arg(&output_masked);
    cmd2.assert().success();
    
    let content_masked = fs::read_to_string(&output_masked)?;
    assert!(content_masked.contains("acgtACGT"));
    
    // No mask
    let mut cmd3 = Command::cargo_bin("pgr")?;
    cmd3.arg("2bit")
        .arg("tofa")
        .arg(&twobit)
        .arg("--no-mask")
        .arg("-o")
        .arg(&output_unmasked);
    cmd3.assert().success();
    
    let content_unmasked = fs::read_to_string(&output_unmasked)?;
    assert!(content_unmasked.contains("ACGTACGT"));

    Ok(())
}
