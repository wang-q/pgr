use assert_cmd::cargo::cargo_bin_cmd;
use tempfile::NamedTempFile;

#[test]
fn command_fq_to_fa() -> anyhow::Result<()> {
    let mut cmd = cargo_bin_cmd!("pgr");
    let input = "@SEQ_ID\nGATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT\n+\n!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65\n";

    let mut file = NamedTempFile::new()?;
    use std::io::Write;
    file.write_all(input.as_bytes())?;

    cmd.arg("fq")
        .arg("to-fa")
        .arg(file.path())
        .assert()
        .success();

    Ok(())
}

#[test]
fn command_fq_interleave_coverage_gap() -> anyhow::Result<()> {
    // 1. 1 file (FQ) -> Output FA
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("fq")
        .arg("interleave")
        .arg("tests/fastq/R1.fq.gz")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;
    assert!(stdout.starts_with(">"));
    assert!(stdout.contains("/1\n"));
    assert!(stdout.contains("/2\n"));
    assert!(!stdout.contains("\n+\n"));

    // 2. 2 files (FQ) -> Output FA
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("fq")
        .arg("interleave")
        .arg("tests/fastq/R1.fq.gz")
        .arg("tests/fastq/R2.fq.gz")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;
    assert!(stdout.starts_with(">"));
    assert!(stdout.contains("/1\n"));
    assert!(stdout.contains("/2\n"));
    assert!(!stdout.contains("\n+\n"));

    // 3. 2 files (FA) -> Output FQ
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("fq")
        .arg("interleave")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/ufasta.fa")
        .arg("--fq")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;
    assert!(stdout.starts_with("@"));
    assert!(stdout.contains("\n+\n"));
    // FA -> FQ fills quality with '!'
    assert!(stdout.contains("!"));

    Ok(())
}

#[test]
fn command_fq_to_fa_output() -> anyhow::Result<()> {
    let mut cmd = cargo_bin_cmd!("pgr");
    let input = "@SEQ_ID\nGATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT\n+\n!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65\n";

    let mut file = NamedTempFile::new()?;
    use std::io::Write;
    file.write_all(input.as_bytes())?;

    let output = cmd
        .arg("fq")
        .arg("to-fa")
        .arg(file.path())
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains(">SEQ_ID"));
    assert!(stdout.contains("GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT"));

    Ok(())
}

#[test]
fn command_fq_to_fa_r1() -> anyhow::Result<()> {
    // Basic conversion test
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("fq")
        .arg("to-fa")
        .arg("tests/fastq/R1.fq.gz")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Verify output format
    assert_eq!(stdout.lines().filter(|e| e.starts_with(">")).count(), 25);
    assert_eq!(stdout.lines().filter(|e| e.is_empty()).count(), 0);
    assert_eq!(stdout.lines().filter(|e| *e == "+").count(), 0);
    assert_eq!(stdout.lines().filter(|e| *e == "!").count(), 0);

    // Test file output
    let mut cmd = cargo_bin_cmd!("pgr");
    let temp = tempfile::Builder::new().suffix(".fa").tempfile()?;
    let temp_path = temp.path();

    cmd.arg("fq")
        .arg("to-fa")
        .arg("tests/fastq/R1.fq.gz")
        .arg("-o")
        .arg(temp_path)
        .assert()
        .success();

    // Read and verify output file
    let output = std::fs::read_to_string(temp_path)?;
    assert_eq!(output.lines().filter(|e| e.starts_with(">")).count(), 25);

    Ok(())
}

#[test]
fn command_fq_interleave() -> anyhow::Result<()> {
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("fq")
        .arg("interleave")
        .arg("tests/fastq/R1.fq.gz")
        .arg("tests/fastq/R2.fq.gz")
        .arg("--fq")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Verify output format
    // 25 pairs * 2 reads/pair = 50 reads
    assert_eq!(stdout.lines().filter(|e| e.starts_with("@")).count(), 50);
    // Check if it's FASTQ (has + lines)
    assert!(stdout.contains("\n+\n"));

    Ok(())
}

#[test]
fn command_fq_interleave_fa() -> anyhow::Result<()> {
    // count empty seqs
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("fq")
        .arg("interleave")
        .arg("tests/fasta/ufasta.fa.gz")
        .arg("tests/fasta/ufasta.fa")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().filter(|e| e.is_empty()).count(), 10);

    // count empty seqs (single)
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("fq")
        .arg("interleave")
        .arg("tests/fasta/ufasta.fa")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().filter(|e| e.is_empty()).count(), 5);

    // count empty seqs (single)
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("fq")
        .arg("interleave")
        .arg("tests/fasta/ufasta.fa")
        .arg("--fq")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().filter(|e| e.is_empty()).count(), 10);

    Ok(())
}

#[test]
fn command_fq_interleave_fq_detailed() -> anyhow::Result<()> {
    // fq
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("fq")
        .arg("interleave")
        .arg("--fq")
        .arg("tests/fastq/R1.fq.gz")
        .arg("tests/fastq/R2.fq.gz")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().filter(|e| *e == "!").count(), 0);
    assert_eq!(stdout.lines().filter(|e| *e == "+").count(), 50);
    assert_eq!(stdout.lines().filter(|e| e.ends_with("/1")).count(), 25);
    assert_eq!(stdout.lines().filter(|e| e.ends_with("/2")).count(), 25);

    // fq (single)
    let mut cmd = cargo_bin_cmd!("pgr");
    let output = cmd
        .arg("fq")
        .arg("interleave")
        .arg("--fq")
        .arg("tests/fastq/R1.fq.gz")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().filter(|e| *e == "!").count(), 25);
    assert_eq!(stdout.lines().filter(|e| *e == "+").count(), 50);
    assert_eq!(stdout.lines().filter(|e| e.ends_with("/1")).count(), 25);
    assert_eq!(stdout.lines().filter(|e| e.ends_with("/2")).count(), 25);

    Ok(())
}
