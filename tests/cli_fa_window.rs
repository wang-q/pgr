use assert_cmd::cargo::cargo_bin_cmd;
use predicates::prelude::*;
use tempfile::TempDir;

#[test]
fn test_fa_window_basic() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = cargo_bin_cmd!("pgr");
    let temp_dir = TempDir::new()?;
    let output_file = temp_dir.path().join("output.fa");

    // Create a simple test file
    let input_file = temp_dir.path().join("input.fa");
    std::fs::write(
        &input_file,
        ">seq1\nATGCATGCAT\n>seq2\nGGCCGGCCGG\n"
    )?;

    cmd.arg("fa")
        .arg("window")
        .arg(&input_file)
        .arg("-l")
        .arg("4")
        .arg("-s")
        .arg("2")
        .arg("-o")
        .arg(&output_file);

    cmd.assert().success();

    // Verify output content
    let content = std::fs::read_to_string(&output_file)?;
    // seq1 (10bp): 1-4, 3-6, 5-8, 7-10, 9-10 -> 5 windows
    // seq2 (10bp): 1-4, 3-6, 5-8, 7-10, 9-10 -> 5 windows
    // Total 10 windows
    assert_eq!(content.matches(">seq1").count(), 5);
    assert_eq!(content.matches(">seq2").count(), 5);
    
    // Check first window
    assert!(content.contains(">seq1:1-4"));
    assert!(content.contains("ATGC"));
    // Check overlapping window
    assert!(content.contains(">seq1:3-6"));
    assert!(content.contains("GCAT"));
    // Check last window (partial)
    assert!(content.contains(">seq1:9-10"));
    assert!(content.contains("AT"));

    Ok(())
}

#[test]
fn test_fa_window_skip_n() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = cargo_bin_cmd!("pgr");
    let temp_dir = TempDir::new()?;
    let output_file = temp_dir.path().join("output.fa");

    let input_file = temp_dir.path().join("input.fa");
    std::fs::write(
        &input_file,
        ">seq1\nATGC\nNNNN\nGCAT\n" // ATGC NNNN GCAT
    )?;

    cmd.arg("fa")
        .arg("window")
        .arg(&input_file)
        .arg("-l")
        .arg("4")
        .arg("-s")
        .arg("4")
        .arg("-o")
        .arg(&output_file);

    cmd.assert().success();

    let content = std::fs::read_to_string(&output_file)?;
    // Should contain ATGC (1-4) and GCAT (9-12), but skip NNNN (5-8)
    assert!(content.contains(">seq1:1-4"));
    assert!(content.contains(">seq1:9-12"));
    assert!(!content.contains(">seq1:5-8"));

    Ok(())
}

#[test]
fn test_fa_window_chunk() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = cargo_bin_cmd!("pgr");
    let temp_dir = TempDir::new()?;
    let output_file = temp_dir.path().join("output.fa");

    let input_file = temp_dir.path().join("input.fa");
    // Create 10 records of 10bp
    let mut input_content = String::new();
    for i in 0..10 {
        input_content.push_str(&format!(">seq{}\nAAAAAAAAAA\n", i));
    }
    std::fs::write(&input_file, input_content)?;

    // Window size 10, step 10 -> 1 window per sequence -> 10 output records total
    // Chunk size 3 -> Should produce 4 files: .001 (3), .002 (3), .003 (3), .004 (1)
    cmd.arg("fa")
        .arg("window")
        .arg(&input_file)
        .arg("-l")
        .arg("10")
        .arg("-s")
        .arg("10")
        .arg("--chunk")
        .arg("3")
        .arg("-o")
        .arg(&output_file);

    cmd.assert().success();

    assert!(temp_dir.path().join("output.001.fa").exists());
    assert!(temp_dir.path().join("output.002.fa").exists());
    assert!(temp_dir.path().join("output.003.fa").exists());
    assert!(temp_dir.path().join("output.004.fa").exists());
    assert!(!temp_dir.path().join("output.005.fa").exists());

    let f1 = std::fs::read_to_string(temp_dir.path().join("output.001.fa"))?;
    assert_eq!(f1.matches(">").count(), 3);
    
    let f4 = std::fs::read_to_string(temp_dir.path().join("output.004.fa"))?;
    assert_eq!(f4.matches(">").count(), 1);

    Ok(())
}

#[test]
fn test_fa_window_shuffle_chunk() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = cargo_bin_cmd!("pgr");
    let temp_dir = TempDir::new()?;
    let output_file = temp_dir.path().join("output.fa");

    let input_file = temp_dir.path().join("input.fa");
    // 100 sequences to ensure shuffle is noticeable (though we check logic mostly)
    let mut input_content = String::new();
    for i in 0..100 {
        input_content.push_str(&format!(">seq{}\nAAAAAAAAAA\n", i));
    }
    std::fs::write(&input_file, input_content)?;

    // Chunk 20 -> 5 files
    cmd.arg("fa")
        .arg("window")
        .arg(&input_file)
        .arg("-l")
        .arg("10")
        .arg("-s")
        .arg("10")
        .arg("--chunk")
        .arg("20")
        .arg("--shuffle")
        .arg("--seed")
        .arg("42")
        .arg("-o")
        .arg(&output_file);

    cmd.assert().success();

    for i in 1..=5 {
        let p = temp_dir.path().join(format!("output.{:03}.fa", i));
        assert!(p.exists(), "File {} should exist", p.display());
        let c = std::fs::read_to_string(p)?;
        assert_eq!(c.matches(">").count(), 20);
    }

    Ok(())
}

#[test]
fn test_fa_window_real_file() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = cargo_bin_cmd!("pgr");
    let temp_dir = TempDir::new()?;
    let output_file = temp_dir.path().join("sakai_window.fa");
    
    // Use the real sakai.fa.gz file
    let input_file = "tests/genome/sakai.fa.gz";

    cmd.arg("fa")
        .arg("window")
        .arg(input_file)
        .arg("-l")
        .arg("1000")
        .arg("-s")
        .arg("500")
        .arg("-o")
        .arg(&output_file);

    cmd.assert().success();

    let content = std::fs::read_to_string(&output_file)?;
    assert!(content.len() > 0);
    assert!(content.contains(">"));
    
    // Verify 1-based coordinates in header
    // e.g., >NC_002695:1-1000
    assert!(content.contains(":1-1000"));

    Ok(())
}

#[test]
fn test_fa_window_chunk_stdout_fail() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = cargo_bin_cmd!("pgr");
    let temp_dir = TempDir::new()?;
    let input_file = temp_dir.path().join("input.fa");
    std::fs::write(&input_file, ">seq1\nACGT\n")?;

    cmd.arg("fa")
        .arg("window")
        .arg(&input_file)
        .arg("--chunk")
        .arg("10"); 
        // Default output is stdout, should fail

    cmd.assert().failure()
        .stderr(predicate::str::contains("Cannot use --chunk with stdout output"));

    Ok(())
}
