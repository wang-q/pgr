use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::fs;
use std::process::Command;
use tempfile::TempDir;

#[test]
fn command_fa_gz() -> anyhow::Result<()> {
    let tempdir = TempDir::new()?;
    let tempdir_str = tempdir.path().to_str().unwrap();

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("fa")
        .arg("gz")
        .arg("tests/index/final.contigs.fa")
        .arg("-o")
        .arg(format!("{}/ctg.fa.gz", tempdir_str))
        .assert()
        .success()
        .stdout(predicate::str::is_empty());

    assert!(&tempdir.path().join("ctg.fa.gz").exists());
    assert!(&tempdir.path().join("ctg.fa.gz.gzi").exists());

    tempdir.close()?;
    Ok(())
}

#[test]
fn command_fa_gz_consistency() -> anyhow::Result<()> {
    if which::which("bgzip").is_err() {
        return Ok(());
    }

    let tempdir = TempDir::new()?;
    let tempdir_str = tempdir.path().to_str().unwrap();
    let infile = "tests/index/final.contigs.fa";
    let outfile_base = "cmp.fa";

    // Run pgr fa gz
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("fa")
        .arg("gz")
        .arg(infile)
        .arg("-o")
        .arg(format!("{}/{}.gz", tempdir_str, outfile_base))
        .assert()
        .success();

    let gz_path = tempdir.path().join(format!("{}.gz", outfile_base));
    let gzi_path = tempdir.path().join(format!("{}.gz.gzi", outfile_base));
    let pgr_gzi_path = tempdir.path().join(format!("{}.gz.gzi.pgr", outfile_base));

    // Rename pgr generated index
    std::fs::rename(&gzi_path, &pgr_gzi_path)?;

    // Run bgzip -r to regenerate index
    let status = std::process::Command::new("bgzip")
        .arg("-r")
        .arg(&gz_path)
        .status()?;
    assert!(status.success());

    // Compare files
    let pgr_gzi = std::fs::read(&pgr_gzi_path)?;
    let bgzip_gzi = std::fs::read(&gzi_path)?;

    assert_eq!(
        pgr_gzi, bgzip_gzi,
        "Generated GZI index differs from bgzip output"
    );

    tempdir.close()?;
    Ok(())
}

#[test]
fn command_fa_gz_level() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let infile = "tests/fasta/ufasta.fa";
    let out_fast = temp.path().join("fast.fa.gz");
    let out_best = temp.path().join("best.fa.gz");

    // Level 1
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("fa")
        .arg("gz")
        .arg(infile)
        .arg("-o")
        .arg(out_fast.to_str().unwrap())
        .arg("-l")
        .arg("1")
        .assert()
        .success();

    // Level 9
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("fa")
        .arg("gz")
        .arg(infile)
        .arg("-o")
        .arg(out_best.to_str().unwrap())
        .arg("-l")
        .arg("9")
        .assert()
        .success();

    // Size check: best should be smaller or equal
    let size_fast = fs::metadata(&out_fast)?.len();
    let size_best = fs::metadata(&out_best)?.len();
    println!("Fast size: {}, Best size: {}", size_fast, size_best);
    assert!(size_best <= size_fast);

    Ok(())
}

#[test]
fn command_fa_gz_reindex() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let raw_infile = "tests/fasta/ufasta.fa";
    let bgzf_file = temp.path().join("work.fa.gz");
    let indexfile = temp.path().join("work.fa.gz.gzi");

    // 1. Create a valid BGZF file first
    let mut cmd_compress = Command::cargo_bin("pgr")?;
    cmd_compress
        .arg("fa")
        .arg("gz")
        .arg(raw_infile)
        .arg("-o")
        .arg(bgzf_file.to_str().unwrap())
        .assert()
        .success();

    assert!(bgzf_file.exists());

    // 2. Run reindex on the generated BGZF file
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("fa")
        .arg("gz")
        .arg(bgzf_file.to_str().unwrap())
        .arg("-r")
        .assert()
        .success();

    assert!(indexfile.exists());

    Ok(())
}

#[test]
fn command_fa_gz_reindex_fail_not_bgzf() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let infile = "tests/fasta/ufasta.fa"; // Normal FASTA, not BGZF

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("fa")
        .arg("gz")
        .arg(infile)
        .arg("-r")
        .assert()
        .failure(); // Should fail

    Ok(())
}

#[test]
fn command_range() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("range")
        .arg("tests/index/final.contigs.fa.gz")
        .arg("k81_130")
        .arg("k81_130:11-20")
        .arg("k81_170:304-323")
        .arg("k81_170(-):1-20")
        .arg("k81_158:70001-70020")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 10);
    assert!(stdout.contains(">k81_130\nAGTTTCAACT"));
    assert!(stdout.contains(">k81_130:11-20\nGGTGAATCAA\n"));
    assert!(stdout.contains(">k81_170:304-323\nAGTTAAAAACCTGATTTATT\n"));
    assert!(stdout.contains(">k81_170(-):1-20\nATTAACCTGTTGTAGGTGTT\n"));
    assert!(stdout.contains(">k81_158:70001-70020\nTGGCTATAACCTAATTTTGT\n"));

    Ok(())
}

#[test]
fn command_range_r() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("range")
        .arg("tests/index/final.contigs.fa.gz")
        .arg("-r")
        .arg("tests/index/sample.rg")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 12);
    assert!(stdout.contains(">k81_130:11-20\nGGTGAATCAA\n"));

    Ok(())
}

#[test]
fn command_range_update() -> anyhow::Result<()> {
    let tempdir = TempDir::new()?;
    let tempdir_str = tempdir.path().to_str().unwrap();

    // Copy test file to temp directory
    std::fs::copy(
        "tests/index/final.contigs.fa",
        format!("{}/test.fa", tempdir_str),
    )?;

    // First run, create index
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("fa")
        .arg("range")
        .arg(format!("{}/test.fa", tempdir_str))
        .arg("k81_130:11-20")
        .output()?;

    // Get index file's modification time
    let loc_file = format!("{}/test.fa.loc", tempdir_str);
    let first_modified = std::fs::metadata(&loc_file)?.modified()?;

    // Wait for a second to ensure different timestamps
    std::thread::sleep(std::time::Duration::from_secs(1));

    // Force update index with --update
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("fa")
        .arg("range")
        .arg(format!("{}/test.fa", tempdir_str))
        .arg("k81_130:11-20")
        .arg("--update")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Verify output content
    assert!(stdout.contains(">k81_130:11-20\nGGTGAATCAA\n"));

    // Verify index file was updated
    let second_modified = std::fs::metadata(&loc_file)?.modified()?;
    assert!(second_modified > first_modified);

    tempdir.close()?;
    Ok(())
}
