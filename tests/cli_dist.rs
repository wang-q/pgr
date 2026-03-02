#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;

#[test]
fn command_dist_hv() {
    let (stdout, _) = PgrCmd::new()
        .args(&["dist", "hv", "tests/clust/IBPA.fa", "-k", "7", "-w", "1"])
        .run();

    assert!(stdout.lines().count() >= 1);
    assert!(stdout.contains("tests/clust/IBPA.fa"));
}

#[test]
fn command_dist_hv_pair() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "hv",
            "tests/clust/IBPA.fa",
            "tests/clust/IBPA.fa", // Compare file against itself
        ])
        .run();

    assert!(stdout.contains("tests/clust/IBPA.fa"));
    // Similarity should be 1.0 / Distance 0.0
    // The output format: <file1> <file2> ... <mash_dist> ...
}

#[test]
fn command_dist_vector() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "vector",
            "tests/clust/domain.tsv",
            "--mode",
            "jaccard",
            "--bin",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 100);
    assert!(stdout
        .contains("Acin_baum_1326584_GCF_025854095_1\tAcin_baum_1326584_GCF_025854095_1\t1.0000"));
}

#[test]
fn command_dist_seq() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "seq",
            "tests/clust/IBPA.fa",
            "-k",
            "7",
            "-w",
            "1",
            "--zero",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 100);
    assert!(stdout.contains("IBPA_ECOLI\tIBPA_ECOLI_GA\t0.0669\t0.4556\t0.6260"));
}

#[test]
fn command_dist_seq_sim() {
    let (stdout, _) = PgrCmd::new()
        .args(&[
            "dist",
            "seq",
            "tests/clust/IBPA.fa",
            "-k",
            "7",
            "-w",
            "1",
            "--zero",
            "--sim",
        ])
        .run();

    assert_eq!(stdout.lines().count(), 100);
    // Mash dist 0.0669 -> Sim 1 - 0.0669 = 0.9331
    assert!(stdout.contains("IBPA_ECOLI\tIBPA_ECOLI_GA\t0.9331\t0.4556\t0.6260"));
}

#[test]
fn command_dist_seq_genome() -> anyhow::Result<()> {
    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("dist")
        .arg("seq")
        .arg("tests/genome/sakai.fa.gz")
        .arg("tests/genome/mg1655.fa.gz")
        .arg("-k")
        .arg("21")
        .arg("-w")
        .arg("5")
        .arg("--hasher")
        .arg("mod")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 2);
    assert!(stdout.contains("NC_002695\tNC_000913\t0."));
    assert!(stdout.contains("NC_002128\tNC_000913\t0."));

    Ok(())
}

#[test]
fn command_dist_seq_merge() -> anyhow::Result<()> {
    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd
        .arg("dist")
        .arg("seq")
        .arg("tests/clust/IBPA.fa")
        .arg("-k")
        .arg("7")
        .arg("-w")
        .arg("1")
        .arg("--merge")
        .arg("--hasher")
        .arg("murmur")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 1);
    assert!(stdout.contains("tests/clust/IBPA.fa\ttests/clust/IBPA.fa\t763"));

    Ok(())
}
