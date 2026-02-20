use assert_cmd::cargo::cargo_bin_cmd;
use predicates::prelude::*;
use std::io::Write;
use tempfile::NamedTempFile;

#[test]
fn test_nwk_support() -> anyhow::Result<()> {
    // 1. Create target tree
    let mut target_file = NamedTempFile::new()?;
    writeln!(target_file, "((A,B),(C,D));")?;

    // 2. Create replicate trees
    let mut replicates_file = NamedTempFile::new()?;
    writeln!(replicates_file, "((A,B),(C,D));")?;
    writeln!(replicates_file, "((A,B),(C,D));")?;
    writeln!(replicates_file, "((A,C),(B,D));")?; // different topology

    // 3. Run command (absolute counts)
    let mut cmd = cargo_bin_cmd!("pgr");
    cmd.arg("nwk")
        .arg("support")
        .arg(target_file.path())
        .arg(replicates_file.path());

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("((A,B)2,(C,D)2)3;"));

    // 4. Run command (percent)
    let mut cmd = cargo_bin_cmd!("pgr");
    cmd.arg("nwk")
        .arg("support")
        .arg(target_file.path())
        .arg(replicates_file.path())
        .arg("--percent");

    // 2/3 * 100 = 66
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("((A,B)66,(C,D)66)100;"));

    Ok(())
}
