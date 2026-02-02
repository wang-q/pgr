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
fn command_masked() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("fa").arg("masked").arg("tests/fasta/ufasta.fa").output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("read46:3-4"), "read46");

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

#[test]
fn command_fa_n50() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test.fa");

    // 100, 100, 200, 200, 400
    // Total = 1000
    // N50 = 200 (at 500, we have 400+200 >= 500)
    fs::write(
        &input,
        ">seq1\nN\n>seq2\nN\n>seq3\nNN\n>seq4\nNN\n>seq5\nNNNN\n"
            .replace("N", "N".repeat(100).as_str()),
    )?;

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("fa").arg("n50").arg(&input).output()?;

    let stdout = String::from_utf8(output.stdout)?;
    assert!(stdout.contains("N50\t200\n"));

    Ok(())
}

#[test]
fn command_fa_n50_stats() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test.fa");

    fs::write(
        &input,
        ">seq1\nN\n>seq2\nN\n>seq3\nNN\n>seq4\nNN\n>seq5\nNNNN\n"
            .replace("N", "N".repeat(100).as_str()),
    )?;

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("n50")
        .arg(&input)
        .arg("-S")
        .arg("-A")
        .arg("-C")
        .arg("-H")
        .output()?;

    let stdout = String::from_utf8(output.stdout)?;
    // N50
    assert!(stdout.contains("200\n"));
    // Sum
    assert!(stdout.contains("1000\n"));
    // Avg
    assert!(stdout.contains("200.00\n"));
    // Count
    assert!(stdout.contains("5\n"));

    Ok(())
}


#[test]
fn command_fa_n50_comprehensive() -> anyhow::Result<()> {
    // display header
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("fa").arg("n50").arg("tests/fasta/ufasta.fa").output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 1);
    assert!(stdout.contains("N50\t314"), "line 1");

    // doesn't display header
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("n50")
        .arg("tests/fasta/ufasta.fa")
        .arg("-H")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 1);
    assert!(!stdout.contains("N50\t314"), "line 1");
    assert!(stdout.contains("314"), "line 1");

    // set genome size (NG50)
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("n50")
        .arg("tests/fasta/ufasta.fa")
        .arg("-H")
        .arg("-g")
        .arg("10000")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 1);
    assert!(stdout.contains("297"), "line 1");

    // sum and average of size
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("n50")
        .arg("tests/fasta/ufasta.fa")
        .arg("-H")
        .arg("-S")
        .arg("-A")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 3);
    assert!(stdout.contains("314\n9317\n186.34"), "line 1,2,3");

    // N10, N90, E-size
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("n50")
        .arg("tests/fasta/ufasta.fa")
        .arg("-H")
        .arg("-E")
        .arg("-N")
        .arg("10")
        .arg("-N")
        .arg("90")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 3);
    assert!(stdout.contains("516\n112\n314.70\n"), "line 1,2,3");

    // transposed
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("n50")
        .arg("tests/fasta/ufasta.fa")
        .arg("-E")
        .arg("-N")
        .arg("10")
        .arg("-N")
        .arg("90")
        .arg("-t")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 2);
    assert!(stdout.contains("N10\tN90\tE\n"), "line 1");
    assert!(stdout.contains("516\t112\t314.70\n"), "line 2");

    Ok(())
}
