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
