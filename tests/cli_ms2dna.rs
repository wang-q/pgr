use assert_cmd::Command;
use predicates::prelude::*;

#[test]
fn command_ms2dna_help() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("ms2dna").arg("--help");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("Convert ms output haplotypes"));
    Ok(())
}

#[test]
fn command_ms2dna_basic_stdin() -> anyhow::Result<()> {
    let input = "\
ms 2 1 -r 0 4
//
segsites: 2
positions: 0.25 0.75
01
10
";
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("ms2dna").arg("--seed").arg("42").write_stdin(input);
    
    // With seed 42, we expect deterministic output
    // nsite=4. 
    // Ancestral generation and mutation logic is tested in unit tests,
    // here we verify the CLI pipeline works and produces FASTA format.
    
    let output = cmd.ok()?;
    let stdout = String::from_utf8(output.stdout)?;
    
    // Check for FASTA headers
    assert!(stdout.contains(">S1"));
    assert!(stdout.contains(">S2"));
    // Check sequences are single lines (not wrapped) and length 4
    let lines: Vec<&str> = stdout.lines().collect();
    // >S1
    // SEQ1
    // >S2
    // SEQ2
    assert_eq!(lines.len(), 4);
    assert_eq!(lines[1].len(), 4);
    assert_eq!(lines[3].len(), 4);

    Ok(())
}

#[test]
fn command_ms2dna_custom_gc() -> anyhow::Result<()> {
    let input = "\
ms 1 1 -r 0 100
//
segsites: 0
positions: 
";
    let mut cmd = Command::cargo_bin("pgr")?;
    // High GC content
    cmd.arg("ms2dna").arg("-g").arg("1.0").arg("-s").arg("123").write_stdin(input);
    
    let output = cmd.ok()?;
    let stdout = String::from_utf8(output.stdout)?;
    
    let lines: Vec<&str> = stdout.lines().collect();
    let seq = lines[1];
    // Check that most bases are G or C
    let gc_count = seq.chars().filter(|c| *c == 'G' || *c == 'C').count();
    assert!(gc_count > 90, "Expected high GC content, got {}", gc_count);
    
    Ok(())
}
