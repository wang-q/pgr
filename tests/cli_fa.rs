use assert_cmd::Command;
use std::fs;
use tempfile::TempDir;

#[test]
fn command_fa_size() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test.fa");

    fs::write(&input, ">seq1\nACGT\n>seq2\nACGTACGT\n")?;

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("fa").arg("size").arg(&input).output()?;

    let stdout = String::from_utf8(output.stdout)?;
    assert!(stdout.contains("seq1\t4\n"));
    assert!(stdout.contains("seq2\t8\n"));

    Ok(())
}

#[test]
fn command_fa_size_file() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("fa").arg("size").arg("tests/fasta/ufasta.fa").output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 50);
    assert!(stdout.contains("read0\t359"), "read0");
    assert!(stdout.contains("read1\t106"), "read1");

    let mut sum = 0;
    for line in stdout.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() == 2 {
            sum += fields[1].parse::<i32>()?;
        }
    }
    assert_eq!(sum, 9317, "sum length");

    Ok(())
}
