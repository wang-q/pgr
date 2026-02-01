use assert_cmd::Command;
use predicates::prelude::*;
use tempfile::tempdir;
use std::io::Write;
use std::fs::File;

#[test]
fn command_axt_help() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("axt").arg("--help");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("Axt tools"));
    Ok(())
}

#[test]
fn command_axt_sort_help() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("axt").arg("sort").arg("--help");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("Sort axt files"));
    Ok(())
}

#[test]
fn command_axt_tomaf_help() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("axt").arg("tomaf").arg("--help");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("Convert from axt to maf format"));
    Ok(())
}

#[test]
fn command_axt_sort_basic() -> anyhow::Result<()> {
    let dir = tempdir()?;
    let input_path = dir.path().join("input.axt");
    let output_path = dir.path().join("output.axt");

    let input_content = "\
0 chr1 11 21 chr2 11 21 - 100
ACTG
ACTG

1 chr1 6 16 chr2 31 41 + 50
AAAA
AAAA

2 chr1 31 41 chr2 6 16 + 200
TTTT
TTTT
";
    {
        let mut f = File::create(&input_path)?;
        f.write_all(input_content.as_bytes())?;
    }

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("axt").arg("sort")
       .arg(&input_path)
       .arg(&output_path);
    
    cmd.assert().success();

    let output = std::fs::read_to_string(&output_path)?;
    let lines: Vec<&str> = output.lines().filter(|l| !l.is_empty()).collect();
    
    // Expected order: 1 (start 6), 0 (start 11), 2 (start 31)
    assert!(lines[0].contains("chr1 6 16"));
    assert!(lines[3].contains("chr1 11 21"));
    assert!(lines[6].contains("chr1 31 41"));

    Ok(())
}

#[test]
fn command_axt_sort_by_score() -> anyhow::Result<()> {
    let dir = tempdir()?;
    let input_path = dir.path().join("input.axt");
    let output_path = dir.path().join("output.axt");

    let input_content = "\
0 chr1 11 21 chr2 11 21 - 100
ACTG
ACTG

1 chr1 6 16 chr2 31 41 + 50
AAAA
AAAA

2 chr1 31 41 chr2 6 16 + 200
TTTT
TTTT
";
    {
        let mut f = File::create(&input_path)?;
        f.write_all(input_content.as_bytes())?;
    }

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("axt").arg("sort")
       .arg("--by-score")
       .arg(&input_path)
       .arg(&output_path);
    
    cmd.assert().success();

    let output = std::fs::read_to_string(&output_path)?;
    let lines: Vec<&str> = output.lines().filter(|l| !l.is_empty()).collect();
    
    // Expected order: 2 (score 200), 0 (score 100), 1 (score 50)
    assert!(lines[0].contains("2 chr1 31 41"));
    assert!(lines[3].contains("0 chr1 11 21"));
    assert!(lines[6].contains("1 chr1 6 16"));

    Ok(())
}

#[test]
fn command_axt_tomaf_basic() -> anyhow::Result<()> {
    let dir = tempdir()?;
    let axt_path = dir.path().join("input.axt");
    let t_sizes_path = dir.path().join("t.sizes");
    let q_sizes_path = dir.path().join("q.sizes");
    let output_path = dir.path().join("output.maf");

    let axt_content = "\
0 chr1 11 14 chr2 11 14 - 100
ACTG
ACTG
";
    {
        let mut f = File::create(&axt_path)?;
        f.write_all(axt_content.as_bytes())?;
    }
    {
        let mut f = File::create(&t_sizes_path)?;
        f.write_all(b"chr1 1000\n")?;
    }
    {
        let mut f = File::create(&q_sizes_path)?;
        f.write_all(b"chr2 2000\n")?;
    }

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("axt").arg("tomaf")
       .arg(&axt_path)
       .arg(&t_sizes_path)
       .arg(&q_sizes_path)
       .arg(&output_path);

    cmd.assert().success();

    let output = std::fs::read_to_string(&output_path)?;
    assert!(output.contains("scoring=blastz"));
    assert!(output.contains("s chr1"));
    assert!(output.contains("s chr2"));
    
    // AXT: chr1 11 14 (1-based, inclusive). Length 4.
    // MAF: start 10 (0-based), size 4.
    assert!(output.contains("chr1                         10          4"));
    
    Ok(())
}

#[test]
fn command_axt_tomaf_split() -> anyhow::Result<()> {
    let dir = tempdir()?;
    let axt_path = dir.path().join("input.axt");
    let t_sizes_path = dir.path().join("t.sizes");
    let q_sizes_path = dir.path().join("q.sizes");
    let output_dir = dir.path().join("output_split");

    let axt_content = "\
0 chr1 11 14 chr2 11 14 - 100
ACTG
ACTG

1 chr2 21 24 chr1 21 24 + 100
ACTG
ACTG
";
    {
        let mut f = File::create(&axt_path)?;
        f.write_all(axt_content.as_bytes())?;
    }
    {
        let mut f = File::create(&t_sizes_path)?;
        f.write_all(b"chr1 1000\nchr2 1000\n")?;
    }
    {
        let mut f = File::create(&q_sizes_path)?;
        f.write_all(b"chr1 2000\nchr2 2000\n")?;
    }

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("axt").arg("tomaf")
       .arg("--t-split")
       .arg(&axt_path)
       .arg(&t_sizes_path)
       .arg(&q_sizes_path)
       .arg(&output_dir);

    cmd.assert().success();

    assert!(output_dir.exists());
    assert!(output_dir.join("chr1.maf").exists());
    assert!(output_dir.join("chr2.maf").exists());

    let output_chr1 = std::fs::read_to_string(output_dir.join("chr1.maf"))?;
    assert!(output_chr1.contains("s chr1                         10"));
    assert!(!output_chr1.contains("s chr2                         20")); // Should not contain the second record

    let output_chr2 = std::fs::read_to_string(output_dir.join("chr2.maf"))?;
    assert!(output_chr2.contains("s chr2                         20"));
    assert!(!output_chr2.contains("s chr1                         10")); // Should not contain the first record

    Ok(())
}
