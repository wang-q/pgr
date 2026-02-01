use assert_cmd::Command;
use flate2::read::GzDecoder;
use predicates::prelude::*;
use std::fs;
use std::io::Read;
use std::path::PathBuf;

#[test]
fn test_lav_to_psl() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let input = r#"#:lav
s {
    "/path/target.fa" 1 1000
    "/path/query.fa" 1 500
}
h {
    ">target.fa"
    ">query.fa"
}
a {
    s 100
    l 1 1 10 10 95
}
"#;

    cmd.arg("lav")
        .arg("topsl")
        .write_stdin(input)
        .assert()
        .success()
        .stdout(predicate::str::contains(
            "10\t0\t0\t0\t0\t0\t0\t0\t+\tquery\t500\t0\t10\ttarget\t1000\t0\t10\t1\t10,\t0,\t0,",
        ));

    Ok(())
}

#[test]
fn test_lav_to_psl_ucsc_new_style() -> anyhow::Result<()> {
    run_ucsc_test("newStyleLastz")
}

#[test]
fn test_lav_to_psl_ucsc_old_style() -> anyhow::Result<()> {
    run_ucsc_test("oldStyleBlastz")
}

#[test]
fn test_lav_to_psl_trim_ends_bug() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let input_path = PathBuf::from("tests/lav/trimEndsBug.lav");

    if !input_path.exists() {
        eprintln!(
            "Skipping test trimEndsBug because input file not found at {:?}",
            input_path
        );
        return Ok(());
    }

    cmd.arg("lav")
        .arg("topsl")
        .arg(input_path)
        .assert()
        .success();

    Ok(())
}

#[test]
fn test_lav_to_psl_redmine12502() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let input_path = PathBuf::from("tests/lav/redmine12502.lav.gz");
    let expected_path = PathBuf::from("tests/lav/redmine12502.psl.gz");

    if !input_path.exists() || !expected_path.exists() {
        eprintln!("Skipping test redmine12502 because files not found");
        return Ok(());
    }

    // Decompress expected output
    let file = fs::File::open(&expected_path)?;
    let mut gz = GzDecoder::new(file);
    let mut expected_content = String::new();
    gz.read_to_string(&mut expected_content)?;

    let expected_lines: Vec<String> = expected_content
        .lines()
        .filter(|l| !l.starts_with("#") && !l.trim().is_empty())
        .map(|l| l.to_string())
        .collect();

    let assert = cmd
        .arg("lav")
        .arg("topsl")
        .arg(input_path)
        .assert()
        .success();

    let output = assert.get_output();
    let output_str = std::str::from_utf8(&output.stdout)?;

    let output_lines: Vec<String> = output_str
        .lines()
        .filter(|l| !l.starts_with("#") && !l.trim().is_empty())
        .map(|l| l.to_string())
        .collect();

    assert_eq!(
        output_lines.len(),
        expected_lines.len(),
        "Line count mismatch for redmine12502"
    );

    for (i, (out, exp)) in output_lines.iter().zip(expected_lines.iter()).enumerate() {
        assert_eq!(out, exp, "Mismatch at line {} for redmine12502", i + 1);
    }

    Ok(())
}

fn run_ucsc_test(name: &str) -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let input_path = PathBuf::from("tests/lav").join(format!("{}.lav", name));
    let expected_path = PathBuf::from("tests/lav").join(format!("{}.psl", name));

    if !input_path.exists() {
        eprintln!(
            "Skipping test {} because input file not found at {:?}",
            name, input_path
        );
        return Ok(());
    }

    let expected_content = fs::read_to_string(&expected_path)?;
    let expected_lines: Vec<String> = expected_content
        .lines()
        .filter(|l| !l.starts_with("#") && !l.trim().is_empty())
        .map(|l| l.to_string())
        .collect();

    let assert = cmd
        .arg("lav")
        .arg("topsl")
        .arg(input_path)
        .assert()
        .success();

    let output = assert.get_output();
    let output_str = std::str::from_utf8(&output.stdout)?;

    let output_lines: Vec<String> = output_str
        .lines()
        .filter(|l| !l.starts_with("#") && !l.trim().is_empty())
        .map(|l| l.to_string())
        .collect();

    assert_eq!(
        output_lines.len(),
        expected_lines.len(),
        "Line count mismatch for {}",
        name
    );

    for (i, (out, exp)) in output_lines.iter().zip(expected_lines.iter()).enumerate() {
        assert_eq!(out, exp, "Mismatch at line {} for {}", i + 1, name);
    }

    Ok(())
}
