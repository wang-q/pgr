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
fn command_fa_size_no_ns() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = temp.path().join("test_nons.fa");
    
    // seq1: 12 bases, 4 Ns (ACGT NNNN ACGT) -> 8 bases
    // seq2: 4 bases, 0 Ns -> 4 bases
    fs::write(&input, ">seq1\nACGTNNNNACGT\n>seq2\nACGT\n")?;

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("fa")
        .arg("size")
        .arg(&input)
        .arg("--no-ns")
        .output()?;

    let stdout = String::from_utf8(output.stdout)?;
    assert!(stdout.contains("seq1\t8\n"));
    assert!(stdout.contains("seq2\t4\n"));

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


// faops filter -l 0 -a 10 -z 50 tests/fasta/ufasta.fa stdout
// faops filter -l 0 -a 1 -u <(cat tests/fasta/ufasta.fa tests/fasta/ufasta.fa) stdout
#[test]
fn command_filter() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("filter")
        .arg("tests/fasta/ufasta.fa")
        .arg("-a")
        .arg("10")
        .arg("-z")
        .arg("50")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 12);
    assert!(!stdout.contains(">read0"), "read0");
    assert!(stdout.contains(">read20"), "read20");

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("filter")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/ufasta.fa.gz")
        .arg("--uniq")
        .arg("-a")
        .arg("1")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 90);

    Ok(())
}

#[test]
fn command_filter_fmt() -> anyhow::Result<()> {
    // faops filter -N tests/fasta/filter.fa stdout
    // faops treats '-' as N, which is incorrect
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("filter")
        .arg("tests/fasta/filter.fa")
        .arg("--iupac")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(!stdout.contains(">iupac\nAMRG"), "iupac");
    assert!(stdout.contains(">iupac\nANNG"), "iupac");
    assert!(stdout.contains(">dash\nA-NG"), "dash not changed");

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("filter")
        .arg("tests/fasta/filter.fa")
        .arg("--dash")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(!stdout.contains(">dash\nA-RG"), "dash");
    assert!(stdout.contains(">dash\nARG"), "dash");

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("filter")
        .arg("tests/fasta/filter.fa")
        .arg("--upper")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(!stdout.contains(">upper\nAtcG"), "upper");
    assert!(stdout.contains(">upper\nATCG"), "upper");

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("filter")
        .arg("tests/fasta/filter.fa")
        .arg("--simplify")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(!stdout.contains(">read.1 simplify\nAGGG"), "simplify");
    assert!(stdout.contains(">read\nAGGG"), "simplify");

    Ok(())
}

#[test]
fn command_count() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("fa").arg("count").arg("tests/fasta/ufasta.fa").output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("read45\t0\t0"), "empty");
    assert!(stdout.contains("total\t9317\t2318"), "total");

    Ok(())
}


#[test]
fn command_one() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("one")
        .arg("tests/fasta/ufasta.fa")
        .arg("read12")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 2);
    assert!(stdout.contains("read12\n"), "read12");

    Ok(())
}

#[test]
fn command_order() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("order")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/list.txt")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 4);
    assert!(stdout.contains("read12\n"), "read12");
    assert!(stdout.contains("read0\n"), "read0");

    Ok(())
}

#[test]
fn command_split_name() -> anyhow::Result<()> {
    let tempdir = TempDir::new()?;
    let tempdir_str = tempdir.path().to_str().unwrap();

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("fa")
        .arg("split")
        .arg("name")
        .arg("tests/fasta/ufasta.fa")
        .arg("-o")
        .arg(tempdir_str)
        .assert()
        .success()
        .stdout(predicates::str::is_empty());

    assert!(&tempdir.path().join("read0.fa").is_file());
    assert!(!&tempdir.path().join("000.fa").exists());

    tempdir.close()?;
    Ok(())
}

#[test]
fn command_split_about() -> anyhow::Result<()> {
    let tempdir = TempDir::new()?;
    let tempdir_str = tempdir.path().to_str().unwrap();

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("fa")
        .arg("split")
        .arg("about")
        .arg("tests/fasta/ufasta.fa")
        .arg("-c")
        .arg("2000")
        .arg("-o")
        .arg(tempdir_str)
        .assert()
        .success()
        .stdout(predicates::str::is_empty());

    assert!(!&tempdir.path().join("read0.fa").is_file());
    assert!(&tempdir.path().join("000.fa").exists());
    assert!(&tempdir.path().join("004.fa").exists());
    assert!(!&tempdir.path().join("005.fa").exists());

    tempdir.close()?;
    Ok(())
}

#[test]
fn command_replace() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("replace")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/replace.tsv")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 95);
    assert!(stdout.contains(">359"), "read0");
    assert!(!stdout.contains(">read0"), "read0");

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("replace")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/replace.tsv")
        .arg("--some")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 6);
    assert!(stdout.contains(">359"), "read0");
    assert!(!stdout.contains(">read0"), "read0");

    Ok(())
}

#[test]
fn command_rc() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("fa").arg("rc").arg("tests/fasta/ufasta.fa").output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("GgacTgcggCTagAA"), "read46");

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("rc")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/list.txt")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains(">RC_read12"), "read12");
    assert!(!stdout.contains(">RC_read46"), "read46");
    assert!(!stdout.contains("GgacTgcggCTagAA"), "read46");

    Ok(())
}

#[test]
fn command_dedup() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("dedup")
        .arg("tests/fasta/dedup.fa")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 8);
    assert!(!stdout.contains(">read0 some text"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("dedup")
        .arg("tests/fasta/dedup.fa")
        .arg("--desc")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 10);
    assert!(stdout.contains(">read0 some text"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("dedup")
        .arg("tests/fasta/dedup.fa")
        .arg("--seq")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 6);
    assert!(!stdout.contains(">read1"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("dedup")
        .arg("tests/fasta/dedup.fa")
        .arg("--seq")
        .arg("--case")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 4);
    assert!(!stdout.contains(">read2"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("dedup")
        .arg("tests/fasta/dedup.fa")
        .arg("--seq")
        .arg("--both")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 2);
    assert!(!stdout.contains(">read3"));

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("dedup")
        .arg("tests/fasta/dedup.fa")
        .arg("--seq")
        .arg("--both")
        .arg("--file")
        .arg("stdout")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 7);
    assert!(stdout.contains(">read0"));
    assert!(stdout.contains("read0\tread3"));

    Ok(())
}


#[test]
fn command_mask() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("mask")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/mask.json")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("read0\ntcgtttaacccaaatcaagg"), "read0");
    assert!(stdout.contains("read2\natagcaagct"), "read2");

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("mask")
        .arg("--hard")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/mask.json")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("read0\nNNNNNNNNNNNNNNNNNNNN"), "read0");
    assert!(stdout.contains("read2\nNNNNNNNNNN"), "read2");

    Ok(())
}
