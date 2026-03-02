#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use tempfile::NamedTempFile;

#[test]
fn command_fq_to_fa() {
    let input = "@SEQ_ID\nGATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT\n+\n!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65\n";

    let mut file = NamedTempFile::new().unwrap();
    use std::io::Write;
    file.write_all(input.as_bytes()).unwrap();

    PgrCmd::new()
        .args(&["fq", "to-fa", file.path().to_str().unwrap()])
        .assert()
        .success();
}

#[test]
fn command_fq_interleave_coverage_gap() {
    // 1. 1 file (FQ) -> Output FA
    let (stdout, _) = PgrCmd::new()
        .args(&["fq", "interleave", "tests/fastq/R1.fq.gz"])
        .run();
    assert!(stdout.starts_with(">"));
    assert!(stdout.contains("/1\n"));
    assert!(stdout.contains("/2\n"));
    assert!(!stdout.contains("\n+\n"));

    // 2. 2 files (FQ) -> Output FA
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fq",
            "interleave",
            "tests/fastq/R1.fq.gz",
            "tests/fastq/R2.fq.gz",
        ])
        .run();
    assert!(stdout.starts_with(">"));
    assert!(stdout.contains("/1\n"));
    assert!(stdout.contains("/2\n"));
    assert!(!stdout.contains("\n+\n"));

    // 3. 2 files (FA) -> Output FQ
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fq",
            "interleave",
            "tests/fasta/ufasta.fa",
            "tests/fasta/ufasta.fa",
            "--fq",
        ])
        .run();
    assert!(stdout.starts_with("@"));
    assert!(stdout.contains("\n+\n"));
    // FA -> FQ fills quality with '!'
    assert!(stdout.contains("!"));
}

#[test]
fn command_fq_to_fa_output() {
    let input = "@SEQ_ID\nGATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT\n+\n!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65\n";

    let mut file = NamedTempFile::new().unwrap();
    use std::io::Write;
    file.write_all(input.as_bytes()).unwrap();

    let (stdout, _) = PgrCmd::new()
        .args(&["fq", "to-fa", file.path().to_str().unwrap()])
        .run();

    assert!(stdout.contains(">SEQ_ID"));
    assert!(stdout.contains("GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT"));
}

#[test]
fn command_fq_to_fa_r1() {
    // Basic conversion test
    let (stdout, _) = PgrCmd::new()
        .args(&["fq", "to-fa", "tests/fastq/R1.fq.gz"])
        .run();

    // Verify output format
    assert_eq!(stdout.lines().filter(|e| e.starts_with(">")).count(), 25);
    assert_eq!(stdout.lines().filter(|e| e.is_empty()).count(), 0);
    assert_eq!(stdout.lines().filter(|e| *e == "+").count(), 0);
    assert_eq!(stdout.lines().filter(|e| *e == "!").count(), 0);

    // Test file output
    let temp = tempfile::Builder::new().suffix(".fa").tempfile().unwrap();
    let temp_path = temp.path();

    PgrCmd::new()
        .args(&[
            "fq",
            "to-fa",
            "tests/fastq/R1.fq.gz",
            "-o",
            temp_path.to_str().unwrap(),
        ])
        .assert()
        .success();

    // Read and verify output file
    let output = std::fs::read_to_string(temp_path).unwrap();
    assert_eq!(output.lines().filter(|e| e.starts_with(">")).count(), 25);
}

#[test]
fn command_fq_interleave() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fq",
            "interleave",
            "tests/fastq/R1.fq.gz",
            "tests/fastq/R2.fq.gz",
            "--fq",
        ])
        .run();

    // Verify output format
    // 25 pairs * 2 reads/pair = 50 reads
    assert_eq!(stdout.lines().filter(|e| e.starts_with("@")).count(), 50);
    // Check if it's FASTQ (has + lines)
    assert!(stdout.contains("\n+\n"));
}

#[test]
fn command_fq_interleave_fa() {
    // count empty seqs
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fq",
            "interleave",
            "tests/fasta/ufasta.fa.gz",
            "tests/fasta/ufasta.fa",
        ])
        .run();

    assert_eq!(stdout.lines().filter(|e| e.is_empty()).count(), 10);

    // count empty seqs (single)
    let (stdout, _) = PgrCmd::new()
        .args(&["fq", "interleave", "tests/fasta/ufasta.fa"])
        .run();

    assert_eq!(stdout.lines().filter(|e| e.is_empty()).count(), 5);

    // count empty seqs (single)
    let (stdout, _) = PgrCmd::new()
        .args(&["fq", "interleave", "tests/fasta/ufasta.fa", "--fq"])
        .run();

    assert_eq!(stdout.lines().filter(|e| e.is_empty()).count(), 10);
}

#[test]
fn command_fq_interleave_fq_detailed() {
    // fq
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "fq",
            "interleave",
            "--fq",
            "tests/fastq/R1.fq.gz",
            "tests/fastq/R2.fq.gz",
        ])
        .run();

    assert_eq!(stdout.lines().filter(|e| *e == "!").count(), 0);
    assert_eq!(stdout.lines().filter(|e| *e == "+").count(), 50);
    assert_eq!(stdout.lines().filter(|e| e.ends_with("/1")).count(), 25);
    assert_eq!(stdout.lines().filter(|e| e.ends_with("/2")).count(), 25);

    // fq (single)
    let (stdout, _) = PgrCmd::new()
        .args(&["fq", "interleave", "--fq", "tests/fastq/R1.fq.gz"])
        .run();

    assert_eq!(stdout.lines().filter(|e| *e == "!").count(), 25);
    assert_eq!(stdout.lines().filter(|e| *e == "+").count(), 50);
    assert_eq!(stdout.lines().filter(|e| e.ends_with("/1")).count(), 25);
    assert_eq!(stdout.lines().filter(|e| e.ends_with("/2")).count(), 25);
}
