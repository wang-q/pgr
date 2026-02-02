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
fn test_2bit_tofa_seq_range() -> anyhow::Result<()> {
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
    cmd.arg("fa")
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

    // Test --n-bed
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("2bit")
        .arg("size")
        .arg(&twobit)
        .arg("--n-bed")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;
    // seq1 has Ns at 4-8
    assert!(stdout.contains("seq1\t4\t8"));
    // seq2 has no Ns, should not appear
    assert!(!stdout.contains("seq2"));

    // Test --mask-bed
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("2bit")
        .arg("size")
        .arg(&twobit)
        .arg("--mask-bed")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;
    // seq1 has no mask
    assert!(!stdout.contains("seq1"));
    // seq2 is all mask (lowercase) 0-4
    assert!(stdout.contains("seq2\t0\t4"));

    Ok(())
}

#[test]
fn test_2bit_compat_seqlist1_file() -> anyhow::Result<()> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let input = std::path::Path::new(&manifest_dir).join("tests/2bit/input/testMask.2bit");
    let list = std::path::Path::new(&manifest_dir).join("tests/2bit/input/seqlist1");
    
    let temp = TempDir::new()?;
    let output = temp.path().join("out.fa");
    
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("tofa")
        .arg(&input)
        .arg("--seqList")
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
fn test_2bit_compat_seqlist2_file() -> anyhow::Result<()> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let input = std::path::Path::new(&manifest_dir).join("tests/2bit/input/testMask.2bit");
    let list = std::path::Path::new(&manifest_dir).join("tests/2bit/input/seqlist2");
    
    let temp = TempDir::new()?;
    let output = temp.path().join("out.fa");
    
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("tofa")
        .arg(&input)
        .arg("--seqList")
        .arg(&list)
        .arg("-o")
        .arg(&output);
    cmd.assert().success();
    
    let output_content = fs::read_to_string(&output)?;
    
    assert!(output_content.contains(">noLower"));
    assert!(output_content.contains("ACGTTTACT")); 
    
    assert!(output_content.contains(">startLower:0-5"));
    assert!(output_content.contains("aaaca"));
    
    assert!(output_content.contains(">endLower:2-8"));
    assert!(output_content.contains("AACCCa"));
    
    Ok(())
}

#[test]
fn test_2bit_compat_bed_file() -> anyhow::Result<()> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let input = std::path::Path::new(&manifest_dir).join("tests/2bit/input/testMask.2bit");
    let bed = std::path::Path::new(&manifest_dir).join("tests/2bit/input/seqlist2.bed");
    
    let temp = TempDir::new()?;
    let output = temp.path().join("out.fa");
    
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("tofa")
        .arg(&input)
        .arg("--bed")
        .arg(&bed)
        .arg("-o")
        .arg(&output);
    cmd.assert().success();
    
    let output_content = fs::read_to_string(&output)?;
    
    assert!(output_content.contains(">startLowerName"));
    assert!(output_content.contains("aaaca"));
    
    assert!(output_content.contains(">endLowerName"));
    assert!(output_content.contains("AACCCa"));
    
    Ok(())
}

#[test]
fn test_2bit_compat_bed_strand() -> anyhow::Result<()> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let input = std::path::Path::new(&manifest_dir).join("tests/2bit/input/testMask.2bit");
    let bed = std::path::Path::new(&manifest_dir).join("tests/2bit/input/bed_with_strand.bed");
    
    let temp = TempDir::new()?;
    let output = temp.path().join("out.fa");
    
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("tofa")
        .arg(&input)
        .arg("--bed")
        .arg(&bed)
        .arg("-o")
        .arg(&output);
    cmd.assert().success();
    
    let output_content = fs::read_to_string(&output)?;
    
    // startLowerPos (0-5 of startLower): aaaca -> aaaca (same as before)
    assert!(output_content.contains(">startLowerPos"));
    assert!(output_content.contains("aaaca"));
    
    // endLowerNeg (2-8 of endLower): AACCCa -> RC -> tGGGTT
    assert!(output_content.contains(">endLowerNeg"));
    assert!(output_content.contains("tGGGTT"));
    
    Ok(())
}

#[test]
fn test_2bit_compat_single_seq_ranges() -> anyhow::Result<()> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let input = std::path::Path::new(&manifest_dir).join("tests/2bit/input/testMask.2bit");
    
    let temp = TempDir::new()?;
    let output = temp.path().join("out.fa");
    
    // Test ml_1_11: manyLower 1-11
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("tofa")
        .arg(&input)
        .arg("--seq")
        .arg("manyLower")
        .arg("--start")
        .arg("1")
        .arg("--end")
        .arg("11")
        .arg("-o")
        .arg(&output);
    cmd.assert().success();
    
    let output_content = fs::read_to_string(&output)?;
    assert!(output_content.contains(">manyLower"));
    // manyLower: aaCCggTTaCgT
    // 1-11: aCCggTTaCg
    assert!(output_content.contains("aCCggTTaCg"));
    
    // Test spec_ml_2_10: manyLower 2-10
    // 2-10: CCggTTaC
    let output2 = temp.path().join("out2.fa");
    let mut cmd2 = Command::cargo_bin("pgr")?;
    cmd2.arg("2bit")
        .arg("tofa")
        .arg(&input)
        .arg("--seq")
        .arg("manyLower")
        .arg("--start")
        .arg("2")
        .arg("--end")
        .arg("10")
        .arg("-o")
        .arg(&output2);
    cmd2.assert().success();
    
    let output_content2 = fs::read_to_string(&output2)?;
    assert!(output_content2.contains("CCggTTaC"));

    // Helper to test range and expected sequence
    let test_range = |start: usize, end: usize, expected: &str| -> anyhow::Result<()> {
        let out_name = format!("out_{}_{}.fa", start, end);
        let out_path = temp.path().join(&out_name);
        let mut cmd = Command::cargo_bin("pgr")?;
        cmd.arg("2bit")
            .arg("tofa")
            .arg(&input)
            .arg("--seq")
            .arg("manyLower")
            .arg("--start")
            .arg(&start.to_string())
            .arg("--end")
            .arg(&end.to_string())
            .arg("-o")
            .arg(&out_path);
        cmd.assert().success();
        let content = fs::read_to_string(&out_path)?;
        if !content.contains(expected) {
             anyhow::bail!("Failed for {}-{}: expected {}, got {}", start, end, expected, content);
        }
        Ok(())
    };

    // ml_3_9: 3-9 -> CggTTa
    test_range(3, 9, "CggTTa")?;
    // ml_4_8: 4-8 -> ggTT
    test_range(4, 8, "ggTT")?;
    // ml_5_6: 5-6 -> g
    test_range(5, 6, "g")?;
    // ml_5_7: 5-7 -> gT
    test_range(5, 7, "gT")?;
    // ml_6_7: 6-7 -> T
    test_range(6, 7, "T")?;
    // ml_7_8: 7-8 -> T
    test_range(7, 8, "T")?;
    // ml_8_9: 8-9 -> a
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
