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
fn command_fq_to_fa_output_same_as_input_rejected() {
    // Regression: `to-fa` used to open the output writer (truncate) before
    // reading the input, so `-o` pointing at the input silently destroyed it.
    let temp = tempfile::Builder::new().suffix(".fq").tempfile().unwrap();
    let path = temp.path().to_str().unwrap().to_string();
    std::fs::write(&path, "@SEQ\nACGT\n+\n!!!!\n").unwrap();

    PgrCmd::new()
        .args(&["fq", "to-fa", &path, "-o", &path])
        .assert()
        .failure()
        .stderr(predicates::str::contains("is also an input file"));

    // The input must be left intact.
    assert_eq!(
        std::fs::read_to_string(&path).unwrap(),
        "@SEQ\nACGT\n+\n!!!!\n"
    );
}

#[test]
fn command_fq_interleave_output_same_as_input_rejected() {
    // Regression: `interleave` used to open the output writer (truncate)
    // before reading the inputs, so `-o` pointing at an input destroyed it.
    let temp = tempfile::Builder::new().suffix(".fq").tempfile().unwrap();
    let path = temp.path().to_str().unwrap().to_string();
    std::fs::write(&path, "@SEQ\nACGT\n+\n!!!!\n").unwrap();

    PgrCmd::new()
        .args(&["fq", "interleave", &path, "-o", &path])
        .assert()
        .failure()
        .stderr(predicates::str::contains("is also an input file"));

    // The input must be left intact.
    assert_eq!(
        std::fs::read_to_string(&path).unwrap(),
        "@SEQ\nACGT\n+\n!!!!\n"
    );
}

#[test]
fn command_fq_interleave_mismatched_read_counts_rejected() {
    // Regression: two-file interleave used to silently truncate to the shorter
    // file; it must now error instead of dropping reads.
    let r1 = tempfile::Builder::new().suffix(".fq").tempfile().unwrap();
    let r2 = tempfile::Builder::new().suffix(".fq").tempfile().unwrap();
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fq");
    std::fs::write(r1.path(), "@R1\nACGT\n+\n!!!!\n").unwrap();
    std::fs::write(r2.path(), "@R2\nAC\n+\n!!\n@R2b\nACGT\n+\n!!!!\n").unwrap();

    PgrCmd::new()
        .args(&[
            "fq",
            "interleave",
            r1.path().to_str().unwrap(),
            r2.path().to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicates::str::contains(
            "paired files have different numbers of reads",
        ));
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

#[test]
fn command_fq_to_fa_ucsc_real_reads() {
    let (stdout, _) = PgrCmd::new()
        .args(&["fq", "to-fa", "tests/fastq/goodHg19.fastq"])
        .run();

    assert!(stdout.starts_with(">SOLEXA-1GA-2_0048_FC6242L:2:1:4066:1152#0/1"));
    assert!(stdout.contains("TGAATAGCTGGAGGAATGCAGACCTCTG"));
}

#[test]
fn command_fq_to_fa_ucsc_malformed_no_panic() {
    // encodeValidate badDdf/bad.fastq: sequence record missing its description
    // line. Must fail with a friendly error, not panic.
    PgrCmd::new()
        .args(&["fq", "to-fa", "tests/fastq/bad.fastq"])
        .assert()
        .failure()
        .stderr(predicates::str::contains("invalid description prefix"));
}
