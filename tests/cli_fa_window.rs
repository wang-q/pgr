#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use tempfile::TempDir;

#[test]
fn test_fa_window_basic() {
    let temp_dir = TempDir::new().unwrap();
    let output_file = temp_dir.path().join("output.fa");

    // Create a simple test file
    let input_file = temp_dir.path().join("input.fa");
    std::fs::write(&input_file, ">seq1\nATGCATGCAT\n>seq2\nGGCCGGCCGG\n").unwrap();

    PgrCmd::new()
        .args(&[
            "fa",
            "window",
            input_file.to_str().unwrap(),
            "-l",
            "4",
            "-s",
            "2",
            "-o",
            output_file.to_str().unwrap(),
        ])
        .assert()
        .success();

    // Verify output content
    let content = std::fs::read_to_string(&output_file).unwrap();
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
}

#[test]
fn test_fa_window_skip_n() {
    let temp_dir = TempDir::new().unwrap();
    let output_file = temp_dir.path().join("output.fa");

    let input_file = temp_dir.path().join("input.fa");
    std::fs::write(
        &input_file,
        ">seq1\nATGC\nNNNN\nGCAT\n", // ATGC NNNN GCAT
    )
    .unwrap();

    PgrCmd::new()
        .args(&[
            "fa",
            "window",
            input_file.to_str().unwrap(),
            "-l",
            "4",
            "-s",
            "4",
            "-o",
            output_file.to_str().unwrap(),
        ])
        .assert()
        .success();

    let content = std::fs::read_to_string(&output_file).unwrap();
    // Should contain ATGC (1-4) and GCAT (9-12), but skip NNNN (5-8)
    assert!(content.contains(">seq1:1-4"));
    assert!(content.contains(">seq1:9-12"));
    assert!(!content.contains(">seq1:5-8"));
}

#[test]
fn test_fa_window_chunk() {
    let temp_dir = TempDir::new().unwrap();
    let output_file = temp_dir.path().join("output.fa");

    let input_file = temp_dir.path().join("input.fa");
    // Create 10 records of 10bp
    let mut input_content = String::new();
    for i in 0..10 {
        input_content.push_str(&format!(">seq{}\nAAAAAAAAAA\n", i));
    }
    std::fs::write(&input_file, input_content).unwrap();

    // Window size 10, step 10 -> 1 window per sequence -> 10 output records total
    // Chunk size 3 -> Should produce 4 files: .001 (3), .002 (3), .003 (3), .004 (1)
    PgrCmd::new()
        .args(&[
            "fa",
            "window",
            input_file.to_str().unwrap(),
            "-l",
            "10",
            "-s",
            "10",
            "--chunk",
            "3",
            "-o",
            output_file.to_str().unwrap(),
        ])
        .assert()
        .success();

    assert!(temp_dir.path().join("output.001.fa").exists());
    assert!(temp_dir.path().join("output.002.fa").exists());
    assert!(temp_dir.path().join("output.003.fa").exists());
    assert!(temp_dir.path().join("output.004.fa").exists());
    assert!(!temp_dir.path().join("output.005.fa").exists());

    let f1 = std::fs::read_to_string(temp_dir.path().join("output.001.fa")).unwrap();
    assert_eq!(f1.matches(">").count(), 3);

    let f4 = std::fs::read_to_string(temp_dir.path().join("output.004.fa")).unwrap();
    assert_eq!(f4.matches(">").count(), 1);
}

#[test]
fn test_fa_window_shuffle_chunk() {
    let temp_dir = TempDir::new().unwrap();
    let output_file = temp_dir.path().join("output.fa");

    let input_file = temp_dir.path().join("input.fa");
    // 100 sequences to ensure shuffle is noticeable (though we check logic mostly)
    let mut input_content = String::new();
    for i in 0..100 {
        input_content.push_str(&format!(">seq{}\nAAAAAAAAAA\n", i));
    }
    std::fs::write(&input_file, input_content).unwrap();

    // Chunk 20 -> 5 files
    PgrCmd::new()
        .args(&[
            "fa",
            "window",
            input_file.to_str().unwrap(),
            "-l",
            "10",
            "-s",
            "10",
            "--chunk",
            "20",
            "--shuffle",
            "--seed",
            "42",
            "-o",
            output_file.to_str().unwrap(),
        ])
        .assert()
        .success();

    for i in 1..=5 {
        let p = temp_dir.path().join(format!("output.{:03}.fa", i));
        assert!(p.exists(), "File {} should exist", p.display());
        let c = std::fs::read_to_string(p).unwrap();
        assert_eq!(c.matches(">").count(), 20);
    }
}

#[test]
fn test_fa_window_real_file() {
    let temp_dir = TempDir::new().unwrap();
    let output_file = temp_dir.path().join("sakai_window.fa");

    // Use the real sakai.fa.gz file
    let input_file = "tests/genome/sakai.fa.gz";

    PgrCmd::new()
        .args(&[
            "fa",
            "window",
            input_file,
            "-l",
            "1000",
            "-s",
            "500",
            "-o",
            output_file.to_str().unwrap(),
        ])
        .assert()
        .success();

    let content = std::fs::read_to_string(&output_file).unwrap();
    assert!(content.len() > 0);
    assert!(content.contains(">"));

    // Verify 1-based coordinates in header
    // e.g., >NC_002695:1-1000
    assert!(content.contains(":1-1000"));
}

#[test]
fn test_fa_window_chunk_stdout_fail() {
    let temp_dir = TempDir::new().unwrap();
    let input_file = temp_dir.path().join("input.fa");
    std::fs::write(&input_file, ">seq1\nACGT\n").unwrap();

    let (_, stderr) = PgrCmd::new()
        .args(&[
            "fa",
            "window",
            input_file.to_str().unwrap(),
            "--chunk",
            "10",
        ])
        .run_fail();

    assert!(stderr.contains("Cannot use --chunk with stdout output"));
}
