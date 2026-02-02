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

#[test]
fn command_fa_size_gz() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("size")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/ufasta.fa.gz")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 100);
    assert!(stdout.contains("read0\t359"), "read0");
    assert!(stdout.contains("read1\t106"), "read1");

    Ok(())
}

#[test]
fn command_fa_some() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test.fa");
    let list = temp.path().join("list.txt");
    let output = temp.path().join("out.fa");

    fs::write(&input, ">seq1\nACGT\n>seq2\nACGTACGT\n>seq3\nTTTT\n")?;
    fs::write(&list, "seq1\nseq3\n")?;

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("fa")
        .arg("some")
        .arg(&input)
        .arg(&list)
        .arg("-o")
        .arg(&output);
    cmd.assert().success();

    let content = fs::read_to_string(&output)?;
    assert!(content.contains(">seq1"));
    assert!(content.contains(">seq3"));
    assert!(!content.contains(">seq2"));

    Ok(())
}

#[test]
fn command_fa_some_invert() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test.fa");
    let list = temp.path().join("list.txt");
    let output = temp.path().join("out.fa");

    fs::write(&input, ">seq1\nACGT\n>seq2\nACGTACGT\n>seq3\nTTTT\n")?;
    fs::write(&list, "seq1\nseq3\n")?;

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("fa")
        .arg("some")
        .arg(&input)
        .arg(&list)
        .arg("--invert")
        .arg("-o")
        .arg(&output);
    cmd.assert().success();

    let content = fs::read_to_string(&output)?;
    assert!(!content.contains(">seq1"));
    assert!(!content.contains(">seq3"));
    assert!(content.contains(">seq2"));

    Ok(())
}
