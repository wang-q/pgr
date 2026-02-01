use assert_cmd::Command;
use std::path::PathBuf;

fn get_input_path(filename: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("tests/psl/swap/input");
    path.push(filename);
    path
}

fn get_expected_path(filename: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("tests/psl/swap/expected");
    path.push(filename);
    path
}

#[test]
fn test_mrna() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("psl")
        .arg("swap")
        .arg(get_input_path("mrna.psl"))
        .arg("-o")
        .arg("stdout")
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();
    let expected = std::fs::read_to_string(get_expected_path("mrnaTest.psl"))?;

    assert_eq!(stdout.replace("\r\n", "\n"), expected.replace("\r\n", "\n"));
    Ok(())
}

#[test]
fn test_trans() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("psl")
        .arg("swap")
        .arg(get_input_path("trans.psl"))
        .arg("-o")
        .arg("stdout")
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();
    let expected = std::fs::read_to_string(get_expected_path("transTest.psl"))?;

    assert_eq!(stdout.replace("\r\n", "\n"), expected.replace("\r\n", "\n"));
    Ok(())
}

#[test]
fn test_mrna_no_rc() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("psl")
        .arg("swap")
        .arg("--no-rc")
        .arg(get_input_path("mrna.psl"))
        .arg("-o")
        .arg("stdout")
        .output()
        .unwrap();

    let stdout = String::from_utf8(output.stdout).unwrap();
    let expected = std::fs::read_to_string(get_expected_path("mrnaNoRcTest.psl"))?;

    assert_eq!(stdout.replace("\r\n", "\n"), expected.replace("\r\n", "\n"));
    Ok(())
}

