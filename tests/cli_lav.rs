use assert_cmd::cargo::cargo_bin_cmd;
use predicates::prelude::*;
use std::fs;
use std::path::PathBuf;

#[test]
fn test_lav_to_psl() -> anyhow::Result<()> {
    let mut cmd = cargo_bin_cmd!("pgr");
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
        .arg("to-psl")
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
    let mut cmd = cargo_bin_cmd!("pgr");
    let input_path = PathBuf::from("tests/lav/trimEndsBug.lav");

    if !input_path.exists() {
        eprintln!(
            "Skipping test trimEndsBug because input file not found at {:?}",
            input_path
        );
        return Ok(());
    }

    cmd.arg("lav")
        .arg("to-psl")
        .arg(input_path)
        .assert()
        .success();

    Ok(())
}

#[test]
fn test_lav_to_psl_lastz_pgr() -> anyhow::Result<()> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let tests_pgr = std::path::Path::new(&manifest_dir).join("tests/pgr");
    let input_path = tests_pgr.join("lastz.lav");
    let expected_path = tests_pgr.join("lastz.psl");

    if !input_path.exists() || !expected_path.exists() {
        eprintln!("Skipping test lastz_pgr because files not found");
        return Ok(());
    }

    let mut cmd = cargo_bin_cmd!("pgr");
    let assert = cmd
        .arg("lav")
        .arg("to-psl")
        .arg(&input_path)
        .assert()
        .success();

    let output = assert.get_output();
    let output_str = std::str::from_utf8(&output.stdout)?;

    let expected_content = fs::read_to_string(&expected_path)?;

    let output_lines: Vec<String> = output_str
        .lines()
        .filter(|l| !l.starts_with("#") && !l.trim().is_empty())
        .map(|l| l.to_string())
        .collect();

    let expected_lines: Vec<String> = expected_content
        .lines()
        .filter(|l| !l.starts_with("#") && !l.trim().is_empty())
        .map(|l| l.to_string())
        .collect();

    assert_eq!(
        output_lines.len(),
        expected_lines.len(),
        "Line count mismatch for lastz.psl"
    );

    for (i, (out, exp)) in output_lines.iter().zip(expected_lines.iter()).enumerate() {
        assert_eq!(out, exp, "Mismatch at line {} for lastz.psl", i + 1);
    }

    Ok(())
}

fn run_ucsc_test(name: &str) -> anyhow::Result<()> {
    let mut cmd = cargo_bin_cmd!("pgr");
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
        .arg("to-psl")
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
