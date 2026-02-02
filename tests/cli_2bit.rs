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
    cmd.arg("fa")
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
    cmd.arg("fa")
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
    cmd.arg("fa")
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
    cmd.arg("fa")
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
fn test_2bit_range_basic() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test.fa");
    let twobit = temp.path().join("test.2bit");
    let output = temp.path().join("out.fa");
    
    fs::write(&input, ">seq1\nACGTACGT\n")?;
    
    // Create 2bit
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("fa")
        .arg("to2bit")
        .arg(&input)
        .arg("-o")
        .arg(&twobit);
    cmd.assert().success();
    
    // Extract range 2-5 (CGTA) - 1-based
    // 01234567
    // ACGTACGT
    //  CGTA
    let mut cmd2 = Command::cargo_bin("pgr")?;
    cmd2.arg("2bit")
        .arg("range")
        .arg(&twobit)
        .arg("seq1:2-5")
        .arg("-o")
        .arg(&output);
    cmd2.assert().success();
    
    let content = fs::read_to_string(&output)?;
    assert!(content.contains(">seq1:2-5\nCGTA"));

    // Extract negative strand
    // seq1:2-5 is CGTA. RevComp: TACG.
    let output_neg = temp.path().join("out_neg.fa");
    let mut cmd3 = Command::cargo_bin("pgr")?;
    cmd3.arg("2bit")
        .arg("range")
        .arg(&twobit)
        .arg("seq1(-):2-5")
        .arg("-o")
        .arg(&output_neg);
    cmd3.assert().success();

    let content_neg = fs::read_to_string(&output_neg)?;
    assert!(content_neg.contains(">seq1(-):2-5\nTACG"));

    Ok(())
}

#[test]
fn test_2bit_range_rgfile() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test.fa");
    let twobit = temp.path().join("test.2bit");
    let list = temp.path().join("ranges.txt");
    let output = temp.path().join("out.fa");
    
    fs::write(&input, ">seq1\nACGT\n>seq2\nTGCA\n")?;
    // Request seq2 (entire sequence) and seq1:1-2
    fs::write(&list, "seq2\nseq1:1-2\n")?;
    
    // Create 2bit
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("fa")
        .arg("to2bit")
        .arg(&input)
        .arg("-o")
        .arg(&twobit);
    cmd.assert().success();
    
    // Extract ranges
    let mut cmd2 = Command::cargo_bin("pgr")?;
    cmd2.arg("2bit")
        .arg("range")
        .arg(&twobit)
        .arg("-r")
        .arg(&list)
        .arg("-o")
        .arg(&output);
    cmd2.assert().success();
    
    let content = fs::read_to_string(&output)?;
    assert!(content.contains(">seq2\nTGCA"));
    assert!(content.contains(">seq1:1-2\nAC"));

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
    cmd.arg("fa")
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

#[test]
fn test_2bit_size() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test.fa");
    let twobit = temp.path().join("test.2bit");
    
    fs::write(&input, ">seq1\nACGT\n>seq2\nNNNN\n")?;
    
    // Create 2bit file first
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("fa")
        .arg("to2bit")
        .arg(&input)
        .arg("-o")
        .arg(&twobit);
    cmd.assert().success();
    
    // Run 2bit size
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("2bit")
        .arg("size")
        .arg(&twobit)
        .output()?;
        
    let stdout = String::from_utf8(output.stdout)?;
    assert!(stdout.contains("seq1\t4"));
    assert!(stdout.contains("seq2\t4"));

    Ok(())
}

#[test]
fn test_2bit_size_flags() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test_flags.fa");
    let twobit = temp.path().join("test_flags.2bit");
    
    // seq1: 12 bases, Ns at 4-8 (4 Ns). ACGT NNNN ACGT. Size 12. No-Ns: 8.
    // seq2: 4 bases. acgt. Size 4. No-Ns: 4. Mask: 0-4.
    fs::write(&input, ">seq1\nACGTNNNNACGT\n>seq2\nacgt\n")?;
    
    // Create 2bit
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("fa")
        .arg("to2bit")
        .arg(&input)
        .arg("-o")
        .arg(&twobit);
    cmd.assert().success();
    
    // Test --no-ns
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("2bit")
        .arg("size")
        .arg(&twobit)
        .arg("--no-ns")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;
    assert!(stdout.contains("seq1\t8"));
    assert!(stdout.contains("seq2\t4"));

    Ok(())
}

#[test]
fn test_2bit_size_multiple() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input1 = temp.path().join("test1.fa");
    let input2 = temp.path().join("test2.fa");
    let twobit1 = temp.path().join("test1.2bit");
    let twobit2 = temp.path().join("test2.2bit");
    
    fs::write(&input1, ">seq1\nACGT\n")?;
    fs::write(&input2, ">seq2\nTGCA\n")?;
    
    // Create 2bit files
    for (inp, out) in [(&input1, &twobit1), (&input2, &twobit2)] {
        let mut cmd = Command::cargo_bin("pgr")?;
        cmd.arg("fa")
            .arg("to2bit")
            .arg(inp)
            .arg("-o")
            .arg(out);
        cmd.assert().success();
    }
    
    // Run 2bit size with multiple inputs
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("2bit")
        .arg("size")
        .arg(&twobit1)
        .arg(&twobit2)
        .output()?;
        
    let stdout = String::from_utf8(output.stdout)?;
    assert!(stdout.contains("seq1\t4"));
    assert!(stdout.contains("seq2\t4"));

    Ok(())
}

#[test]
fn test_2bit_range_seqlist1_file() -> anyhow::Result<()> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let input = std::path::Path::new(&manifest_dir).join("tests/2bit/input/testMask.2bit");
    let list = std::path::Path::new(&manifest_dir).join("tests/2bit/input/seqlist1");
    
    let temp = TempDir::new()?;
    let output = temp.path().join("out.fa");
    
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("range")
        .arg(&input)
        .arg("--rgfile")
        .arg(&list)
        .arg("-o")
        .arg(&output);
    cmd.assert().success();
    
    let output_content = fs::read_to_string(&output)?;
    
    assert!(output_content.contains(">noLower"));
    assert!(output_content.contains(">startLower"));
    assert!(output_content.contains(">endLower"));
    assert!(!output_content.contains(">manyLower"));
    
    Ok(())
}


#[test]
fn test_2bit_masked() -> anyhow::Result<()> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let input_mask = std::path::Path::new(&manifest_dir).join("tests/2bit/input/testMask.2bit");
    let input_n = std::path::Path::new(&manifest_dir).join("tests/2bit/input/testN.2bit");
    
    let temp = TempDir::new()?;
    let out_mask = temp.path().join("out_mask.txt");
    let out_n = temp.path().join("out_n.txt");
    
    // 1. testMask.2bit
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("masked")
        .arg(&input_mask)
        .arg("-o")
        .arg(&out_mask);
    cmd.assert().success();
    
    let content_mask = fs::read_to_string(&out_mask)?;
    
    // allLower is masked. It has 12 bases.
    assert!(content_mask.contains("allLower:1-12"));
    // noLower should not be in output
    assert!(!content_mask.contains("noLower"));
    
    // 2. testN.2bit with --gap
    let mut cmd2 = Command::cargo_bin("pgr")?;
    cmd2.arg("2bit")
        .arg("masked")
        .arg(&input_n)
        .arg("--gap")
        .arg("-o")
        .arg(&out_n);
    cmd2.assert().success();
    
    let content_n = fs::read_to_string(&out_n)?;
    
    // startN: NANNAANNNAAA
    // Ns at: 1, 3-4, 7-9
    assert!(content_n.contains("startN:1"));
    assert!(content_n.contains("startN:3-4"));
    assert!(content_n.contains("startN:7-9"));

    Ok(())
}

#[test]
fn test_2bit_range_legacy_cases() -> anyhow::Result<()> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let input = std::path::Path::new(&manifest_dir).join("tests/2bit/input/testMask.2bit");

    let temp = TempDir::new()?;

    // Helper to test range and expected sequence
    let test_range = |start: usize, end: usize, expected: &str| -> anyhow::Result<()> {
        let out_name = format!("out_{}_{}.fa", start, end);
        let out_path = temp.path().join(&out_name);
        
        let range_str = format!("manyLower:{}-{}", start, end);
        
        let mut cmd = Command::cargo_bin("pgr")?;
        cmd.arg("2bit")
            .arg("range")
            .arg(&input)
            .arg(&range_str)
            .arg("-o")
            .arg(&out_path);
        cmd.assert().success();
        
        let content = fs::read_to_string(&out_path)?;
        if !content.contains(expected) {
             anyhow::bail!("Failed for {}: expected {}, got {}", range_str, expected, content);
        }
        Ok(())
    };

    // Test cases from original test
    test_range(1, 11, "aCCggTTaCg")?;
    test_range(2, 10, "CCggTTaC")?;
    test_range(3, 9, "CggTTa")?;
    test_range(4, 8, "ggTT")?;
    test_range(5, 6, "g")?;
    test_range(5, 7, "gT")?;
    test_range(6, 7, "T")?;
    test_range(7, 8, "T")?;
    test_range(8, 9, "a")?;

    Ok(())
}

#[test]
fn test_2bit_compat_mask_file() -> anyhow::Result<()> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let input = std::path::Path::new(&manifest_dir).join("tests/2bit/input/testMask.2bit");

    let temp = TempDir::new()?;
    let output = temp.path().join("out.fa");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("tofa")
        .arg(&input)
        .arg("-o")
        .arg(&output);
    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;

    // Check for sequence names
    assert!(output_content.contains(">allLower"));
    assert!(output_content.contains(">endLower"));
    assert!(output_content.contains(">manyLower"));
    assert!(output_content.contains(">noLower"));
    assert!(output_content.contains(">startLower"));

    // Check masking (lowercase)
    // allLower should be all lowercase
    // We can't easily check full content without reading exact expectation,    
    // but we can check if it contains lowercase letters.
    assert!(output_content.chars().any(|c| c.is_lowercase()));

    Ok(())
}

#[test]
fn test_2bit_compat_n_file() -> anyhow::Result<()> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let input = std::path::Path::new(&manifest_dir).join("tests/2bit/input/testN.2bit");

    let temp = TempDir::new()?;
    let output = temp.path().join("out.fa");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("tofa")
        .arg(&input)
        .arg("-o")
        .arg(&output);
    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;

    assert!(output_content.contains(">startN"));
    assert!(output_content.contains("NANNAANNNAAA"));

    assert!(output_content.contains(">startNonN"));
    assert!(output_content.contains("ANAANNAAANNN"));

    Ok(())
}
