use assert_cmd::prelude::*;
use std::process::Command;

#[test]
fn command_dist_hv() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("dist")
        .arg("hv")
        .arg("tests/clust/IBPA.fa")
        .arg("-k")
        .arg("7")
        .arg("-w")
        .arg("1")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.lines().count() >= 1);
    assert!(stdout.contains("tests/clust/IBPA.fa"));

    Ok(())
}

#[test]
fn command_dist_hv_pair() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("dist")
        .arg("hv")
        .arg("tests/clust/IBPA.fa")
        .arg("tests/clust/IBPA.fa") // Compare file against itself
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("tests/clust/IBPA.fa"));
    // Similarity should be 1.0 / Distance 0.0
    // The output format: <file1> <file2> ... <mash_dist> ...

    Ok(())
}

#[test]
fn command_dist_vector() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("dist")
        .arg("vector")
        .arg("tests/clust/domain.tsv")
        .arg("--mode")
        .arg("jaccard")
        .arg("--bin")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 100);
    assert!(stdout
        .contains("Acin_baum_1326584_GCF_025854095_1\tAcin_baum_1326584_GCF_025854095_1\t1.0000"));

    Ok(())
}
