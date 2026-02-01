use assert_cmd::Command;
use std::fs;
use std::path::PathBuf;

fn get_expected_path(filename: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("tests/psl/expected");
    path.push(filename);
    path
}

fn get_input_path(filename: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("tests/psl");
    path.push(filename);
    path
}

fn assert_sorted_eq(stdout: &str, expected: &str) {
    let mut stdout_lines: Vec<&str> = stdout.lines().collect();
    let mut expected_lines: Vec<&str> = expected.lines().collect();

    // Sort from the second line (index 1) onwards to preserve header
    if stdout_lines.len() > 1 {
        stdout_lines[1..].sort();
    }
    if expected_lines.len() > 1 {
        expected_lines[1..].sort();
    }

    assert_eq!(stdout_lines, expected_lines);
}

#[test]
fn untrans_align_test() {
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("psl")
        .arg("stats")
        .arg(get_input_path("stats_basic.psl"))
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();
    let expected = fs::read_to_string(get_expected_path("untransAlignTest.stats")).unwrap();

    assert_sorted_eq(&stdout, &expected);
}

#[test]
fn untrans_query_test() {
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("psl")
        .arg("stats")
        .arg("--query-stats")
        .arg(get_input_path("stats_basic.psl"))
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();
    let expected = fs::read_to_string(get_expected_path("untransQueryTest.stats")).unwrap();

    assert_sorted_eq(&stdout, &expected);
}

#[test]
fn untrans_overall_test() {
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("psl")
        .arg("stats")
        .arg("--overall-stats")
        .arg(get_input_path("stats_basic.psl"))
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();
    let expected = fs::read_to_string(get_expected_path("untransOverallTest.stats")).unwrap();

    assert_sorted_eq(&stdout, &expected);
}

#[test]
fn untrans_align_unaln_test() {
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("psl")
        .arg("stats")
        .arg("--queries")
        .arg(get_input_path("basic.queries"))
        .arg(get_input_path("stats_basic.psl"))
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();
    let expected = fs::read_to_string(get_expected_path("untransAlignUnalnTest.stats")).unwrap();

    assert_sorted_eq(&stdout, &expected);
}

#[test]
fn untrans_query_unaln_test() {
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("psl")
        .arg("stats")
        .arg("--query-stats")
        .arg("--queries")
        .arg(get_input_path("basic.queries"))
        .arg(get_input_path("stats_basic.psl"))
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();
    let expected = fs::read_to_string(get_expected_path("untransQueryUnalnTest.stats")).unwrap();

    assert_sorted_eq(&stdout, &expected);
}

#[test]
fn untrans_overall_unaln_test() {
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("psl")
        .arg("stats")
        .arg("--overall-stats")
        .arg("--queries")
        .arg(get_input_path("basic.queries"))
        .arg(get_input_path("stats_basic.psl"))
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();
    let expected = fs::read_to_string(get_expected_path("untransOverallUnalnTest.stats")).unwrap();

    assert_sorted_eq(&stdout, &expected);
}

#[test]
fn cdna2cdna_align_test() {
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("psl")
        .arg("stats")
        .arg("--queries")
        .arg(get_input_path("cdna2cdna.queries"))
        .arg(get_input_path("cdna2cdna.psl"))
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();
    let expected = fs::read_to_string(get_expected_path("cdna2cdnaAlignTest.stats")).unwrap();

    assert_sorted_eq(&stdout, &expected);
}

#[test]
fn cdna2cdna_query_test() {
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("psl")
        .arg("stats")
        .arg("--query-stats")
        .arg("--queries")
        .arg(get_input_path("cdna2cdna.queries"))
        .arg(get_input_path("cdna2cdna.psl"))
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();
    let expected = fs::read_to_string(get_expected_path("cdna2cdnaQueryTest.stats")).unwrap();

    assert_sorted_eq(&stdout, &expected);
}

#[test]
fn cdna2cdna_overall_test() {
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("psl")
        .arg("stats")
        .arg("--overall-stats")
        .arg("--queries")
        .arg(get_input_path("cdna2cdna.queries"))
        .arg(get_input_path("cdna2cdna.psl"))
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();
    let expected = fs::read_to_string(get_expected_path("cdna2cdnaOverallTest.stats")).unwrap();

    assert_sorted_eq(&stdout, &expected);
}

#[test]
fn tsv_header_test() {
    let mut cmd = Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("psl")
        .arg("stats")
        .arg("--tsv")
        .arg(get_input_path("stats_basic.psl"))
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();
    let expected = fs::read_to_string(get_expected_path("tsvHeaderTest.stats")).unwrap();

    assert_sorted_eq(&stdout, &expected);
}
