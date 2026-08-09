#[macro_use]
#[path = "common/mod.rs"]
mod common;

use std::path::Path;

const HIST_FILE_LEN: u64 = 28 + 32767 * 8;

fn write_fa(path: &Path) {
    std::fs::write(
        path,
        ">chr1\nACGTACGTACGTACGTACGTACGTACGT\n>chr2\nTTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCCAAAAA\n",
    )
    .unwrap();
}

fn write_fq(path: &Path) {
    std::fs::write(
        path,
        "@r1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n@r2\nTTTTGGGGCCCCAAAA\n+\nIIIIIIIIIIIIIIII\n",
    )
    .unwrap();
}

#[test]
fn command_kmer_help() -> anyhow::Result<()> {
    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    let output = cmd.arg("kmer").arg("--help").output()?;
    let stdout = String::from_utf8(output.stdout)?;
    assert!(stdout.contains("Analyzes k-mer counts, profiles"));
    assert!(stdout.contains("table"));
    assert!(stdout.contains("profile"));
    assert!(stdout.contains("hist"));
    Ok(())
}

#[test]
fn command_kmer_table_end_to_end() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let fa = temp.path().join("in.fa");
    let pkt = temp.path().join("t.pkt");
    write_fa(&fa);

    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "kmer",
            "table",
            fa.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            pkt.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("14 unique 8-mers"), "stderr: {stderr}");
    assert!(pkt.exists());

    // The table can be reused: histogram from -t matches the sequence path.
    let h1 = temp.path().join("h1.hist");
    let h2 = temp.path().join("h2.hist");
    common::PgrCmd::new()
        .args(&[
            "kmer",
            "hist",
            fa.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            h1.to_str().unwrap(),
        ])
        .run();
    common::PgrCmd::new()
        .args(&[
            "kmer",
            "hist",
            "-t",
            pkt.to_str().unwrap(),
            "-o",
            h2.to_str().unwrap(),
        ])
        .run();
    assert_eq!(h1.metadata()?.len(), HIST_FILE_LEN);
    assert_eq!(
        std::fs::read(&h1)?,
        std::fs::read(&h2)?,
        "hist from table must match hist from sequences"
    );
    Ok(())
}

#[test]
fn command_kmer_profile_self_and_relative() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let fa = temp.path().join("in.fa");
    let pkt = temp.path().join("t.pkt");
    let self_pkp = temp.path().join("self.pkp");
    let rel_pkp = temp.path().join("rel.pkp");
    write_fa(&fa);

    common::PgrCmd::new()
        .args(&[
            "kmer",
            "table",
            fa.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            pkt.to_str().unwrap(),
        ])
        .run();

    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "kmer",
            "profile",
            fa.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            self_pkp.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("2 profiles"), "stderr: {stderr}");
    assert_eq!(&std::fs::read(&self_pkp)?[0..4], b"PKPP");

    // Relative profile reuses the table; k is read from the table when the
    // command line omits --kmer.
    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "kmer",
            "profile",
            fa.to_str().unwrap(),
            "-t",
            pkt.to_str().unwrap(),
            "-o",
            rel_pkp.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("2 profiles"), "stderr: {stderr}");
    assert_eq!(&std::fs::read(&rel_pkp)?[0..4], b"PKPP");
    Ok(())
}

#[test]
fn command_kmer_reads_fastq_and_stdin() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let fq = temp.path().join("in.fq");
    let pkt = temp.path().join("fq.pkt");
    write_fq(&fq);

    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "kmer",
            "table",
            fq.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            pkt.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("8 unique 8-mers"), "stderr: {stderr}");

    let fa = temp.path().join("in.fa");
    write_fa(&fa);
    let stdin_pkt = temp.path().join("stdin.pkt");
    let input = std::fs::read_to_string(&fa)?;
    let (_, stderr) = common::PgrCmd::new()
        .stdin(input)
        .args(&[
            "kmer",
            "table",
            "stdin",
            "-k",
            "8",
            "-o",
            stdin_pkt.to_str().unwrap(),
        ])
        .run();
    assert!(stderr.contains("14 unique 8-mers"), "stderr: {stderr}");
    Ok(())
}

#[test]
fn command_kmer_argument_validation() -> anyhow::Result<()> {
    let temp = tempfile::TempDir::new()?;
    let fa = temp.path().join("in.fa");
    let pkt = temp.path().join("t.pkt");
    write_fa(&fa);
    common::PgrCmd::new()
        .args(&[
            "kmer",
            "table",
            fa.to_str().unwrap(),
            "-k",
            "8",
            "-o",
            pkt.to_str().unwrap(),
        ])
        .run();

    // --kmer mismatching the table must fail.
    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "kmer",
            "hist",
            "-t",
            pkt.to_str().unwrap(),
            "-k",
            "10",
            "-o",
            temp.path().join("x.hist").to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("does not match table k"),
        "stderr: {stderr}"
    );

    // No --kmer and no --table must fail.
    let (_, stderr) = common::PgrCmd::new()
        .args(&[
            "kmer",
            "hist",
            fa.to_str().unwrap(),
            "-o",
            temp.path().join("y.hist").to_str().unwrap(),
        ])
        .run_fail();
    assert!(stderr.contains("--kmer is required"), "stderr: {stderr}");
    Ok(())
}
