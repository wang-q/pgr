use assert_cmd::Command;
use tempfile::NamedTempFile;

#[test]
fn command_fq_tofa() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let input = "@SEQ_ID\nGATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT\n+\n!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65\n";
    
    let mut file = NamedTempFile::new()?;
    use std::io::Write;
    file.write_all(input.as_bytes())?;
    
    cmd.arg("fq")
        .arg("tofa")
        .arg(file.path())
        .assert()
        .success();

    Ok(())
}

#[test]
fn command_fq_tofa_output() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let input = "@SEQ_ID\nGATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT\n+\n!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65\n";
    
    let mut file = NamedTempFile::new()?;
    use std::io::Write;
    file.write_all(input.as_bytes())?;
    
    let output = cmd.arg("fq")
        .arg("tofa")
        .arg(file.path())
        .output()
        .unwrap();
        
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains(">SEQ_ID"));
    assert!(stdout.contains("GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT"));
    
    Ok(())
}

#[test]
fn command_fq_tofa_r1() -> anyhow::Result<()> {
    // Basic conversion test
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("fq")
        .arg("tofa")
        .arg("tests/fastq/R1.fq.gz")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Verify output format
    assert_eq!(
        stdout
            .lines()
            .into_iter()
            .filter(|e| e.starts_with(">"))
            .count(),
        25
    );
    assert_eq!(stdout.lines().into_iter().filter(|e| *e == "").count(), 0);
    assert_eq!(stdout.lines().into_iter().filter(|e| *e == "+").count(), 0);
    assert_eq!(stdout.lines().into_iter().filter(|e| *e == "!").count(), 0);

    // Test file output
    let mut cmd = Command::cargo_bin("pgr")?;
    let temp = tempfile::Builder::new().suffix(".fa").tempfile()?;
    let temp_path = temp.path();
    
    cmd.arg("fq")
        .arg("tofa")
        .arg("tests/fastq/R1.fq.gz")
        .arg("-o")
        .arg(temp_path)
        .assert()
        .success();

    // Read and verify output file
    let output = std::fs::read_to_string(temp_path)?;
    assert_eq!(
        output
            .lines()
            .into_iter()
            .filter(|e| e.starts_with(">"))
            .count(),
        25
    );

    Ok(())
}
