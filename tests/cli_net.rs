use assert_cmd::cargo::cargo_bin_cmd;
use predicates::prelude::*;
use std::fs;
use std::io::Write;
use tempfile::{NamedTempFile, TempDir};

fn create_2bit(
    dir: &TempDir,
    name: &str,
    content: &str,
) -> Result<std::path::PathBuf, Box<dyn std::error::Error>> {
    let fa_path = dir.path().join(format!("{}.fa", name));
    let bit_path = dir.path().join(format!("{}.2bit", name));
    fs::write(&fa_path, content)?;

    let mut cmd = cargo_bin_cmd!("pgr");
    cmd.arg("fa")
        .arg("to-2bit")
        .arg(&fa_path)
        .arg("-o")
        .arg(&bit_path);
    cmd.assert().success();

    Ok(bit_path)
}

// --- net syntenic tests ---

#[test]
fn test_net_syntenic_basic() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = cargo_bin_cmd!("pgr");

    // Create input net file
    let mut in_file = NamedTempFile::new()?;
    writeln!(in_file, "net chr1 1000")?;
    writeln!(
        in_file,
        " fill 100 200 chr2 + 50 200 id 1 score 100 ali 200"
    )?;
    writeln!(
        in_file,
        " fill 500 100 chr2 + 100 100 id 2 score 50 ali 100"
    )?;

    let out_file = NamedTempFile::new()?;
    let out_path = out_file.path().to_str().unwrap();

    cmd.arg("net")
        .arg("syntenic")
        .arg(in_file.path().to_str().unwrap())
        .arg(out_path)
        .assert()
        .success();

    let output = std::fs::read_to_string(out_path)?;
    println!("Output:\n{}", output);

    // Check output content
    // Fill 1: qDup 100. Type top.
    // Fill 2: qDup 100. Type top.

    assert!(output.contains("fill 100 200 chr2 + 50 200 id 1 score 100 ali 200 qDup 100 type top"));
    assert!(output.contains("fill 500 100 chr2 + 100 100 id 2 score 50 ali 100 qDup 100 type top"));

    Ok(())
}

#[test]
fn test_net_syntenic_nested() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = cargo_bin_cmd!("pgr");

    // Nested structure to test syn/inv/nonSyn
    let mut in_file = NamedTempFile::new()?;
    writeln!(in_file, "net chr1 1000")?;
    writeln!(in_file, " fill 0 500 chr2 + 0 500 id 1 score 100 ali 500")?;
    writeln!(in_file, "  gap 100 100 chr2 + 100 100")?;
    writeln!(
        in_file,
        "   fill 100 100 chr2 + 100 100 id 2 score 50 ali 100"
    )?;

    // Fill 1: 0-500.
    // Gap: 100-200.
    // Fill 2: 100-200.
    // Fill 2 parent is Fill 1.
    // qName same.
    // Intersection: Fill 2 (100-200) vs Fill 1 (0-500). Intersection 100.
    // Strand same (+).
    // Type should be "syn".
    // qOver should be 100.
    // qFar should be 0.

    let out_file = NamedTempFile::new()?;
    let out_path = out_file.path().to_str().unwrap();

    cmd.arg("net")
        .arg("syntenic")
        .arg(in_file.path().to_str().unwrap())
        .arg(out_path)
        .assert()
        .success();

    let output = std::fs::read_to_string(out_path)?;
    println!("Output:\n{}", output);

    assert!(output.contains(
        "fill 100 100 chr2 + 100 100 id 2 score 50 ali 100 qOver 100 qFar 0 qDup 0 type syn"
    ));

    Ok(())
}

// --- net to-axt tests ---

#[test]
fn test_net_to_axt_basic() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = cargo_bin_cmd!("pgr");
    let temp = TempDir::new()?;

    // Create 2bit files
    let t_2bit = create_2bit(&temp, "chrT", ">chrT\nACGTACGTAC")?; // 10 bases
    let q_2bit = create_2bit(&temp, "chrQ", ">chrQ\nACGTACGTAC")?; // 10 bases

    // Create Chain file
    // chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
    // Match 0-10 on T to 0-10 on Q.
    let chain_path = temp.path().join("in.chain");
    let mut chain_file = fs::File::create(&chain_path)?;
    writeln!(chain_file, "chain 100 chrT 10 + 0 10 chrQ 10 + 0 10 1")?;
    writeln!(chain_file, "10")?; // block size 10
    writeln!(chain_file)?;

    // Create Net file
    // net chrT 10
    //  fill 0 10 chrQ + 0 10 id 1 score 100 ali 10
    let net_path = temp.path().join("in.net");
    let mut net_file = fs::File::create(&net_path)?;
    writeln!(net_file, "net chrT 10")?;
    writeln!(net_file, " fill 0 10 chrQ + 0 10 id 1 score 100 ali 10")?;

    let out_path = temp.path().join("out.axt");

    cmd.arg("net")
        .arg("to-axt")
        .arg(&net_path)
        .arg(&chain_path)
        .arg(&t_2bit)
        .arg(&q_2bit)
        .arg(&out_path)
        .assert()
        .success();

    let output = fs::read_to_string(&out_path)?;
    println!("Output:\n{}", output);

    // Check AXT output
    // 0 chrT 1 10 chrQ 1 10 + 100
    // ACGTACGTAC
    // ACGTACGTAC

    assert!(output.contains("0 chrT 1 10 chrQ 1 10 + 955"));
    assert!(output.contains("ACGTACGTAC"));

    Ok(())
}

#[test]
fn test_net_to_axt_reverse() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = cargo_bin_cmd!("pgr");
    let temp = TempDir::new()?;

    // T: ACGT...
    // Q: ACGT...
    // Match T[0-4] to Q[0-4] (on - strand).
    // Q - strand: ACGT... (revcomp of Q +)
    // Q +: ...ACGT
    // Let's make Q+: TTTT ACGT TTTT (len 12)
    // RevComp(Q+): AAAA ACGT AAAA
    // Index 0-4 on - strand is "AAAA".

    let t_2bit = create_2bit(&temp, "chrT", ">chrT\nAAAA")?;
    let q_2bit = create_2bit(&temp, "chrQ", ">chrQ\nTTTTACGTTTTT")?; // 12 bases

    // Chain:
    // Match T:0-4 ("AAAA") to Q-:0-4 ("AAAA").
    // Q-: 0-4 corresponds to Q+: [12-4, 12-0) = [8, 12).
    // Q+: ...TTTT (Wait, last 4 bases of TTTTACGTTTTT is TTTT).
    // RevComp(TTTT) = AAAA. Correct.

    let chain_path = temp.path().join("in.chain");
    let mut chain_file = fs::File::create(&chain_path)?;
    writeln!(chain_file, "chain 100 chrT 4 + 0 4 chrQ 12 - 0 4 1")?;
    writeln!(chain_file, "4")?;
    writeln!(chain_file)?;

    let net_path = temp.path().join("in.net");
    let mut net_file = fs::File::create(&net_path)?;
    writeln!(net_file, "net chrT 4")?;
    writeln!(net_file, " fill 0 4 chrQ - 0 4 id 1 score 100 ali 4")?;

    let out_path = temp.path().join("out.axt");

    cmd.arg("net")
        .arg("to-axt")
        .arg(&net_path)
        .arg(&chain_path)
        .arg(&t_2bit)
        .arg(&q_2bit)
        .arg(&out_path)
        .assert()
        .success();

    let output = fs::read_to_string(&out_path)?;
    println!("Output:\n{}", output);

    // 0 chrT 1 4 chrQ 1 4 - 100
    // AAAA
    // AAAA

    assert!(output.contains("0 chrT 1 4 chrQ 1 4 - 364"));
    assert!(output.contains("AAAA"));
    // Should verify it appears twice (sequence lines)

    Ok(())
}

// --- net split tests ---

#[test]
fn test_net_split() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = cargo_bin_cmd!("pgr");

    // Create input net file with multiple chromosomes
    let mut in_file = NamedTempFile::new()?;
    writeln!(in_file, "net chr1 1000")?;
    writeln!(in_file, " fill 0 100 chrA + 0 100 id 1 score 100 ali 100")?;
    writeln!(in_file, "net chr2 2000")?;
    writeln!(
        in_file,
        " fill 100 200 chrB - 100 200 id 2 score 200 ali 200"
    )?;

    // Create output directory
    let out_dir = tempfile::tempdir()?;
    let out_dir_path = out_dir.path().to_str().unwrap();

    cmd.arg("net")
        .arg("split")
        .arg(in_file.path().to_str().unwrap())
        .arg(out_dir_path)
        .assert()
        .success();

    // Verify output files
    let chr1_path = out_dir.path().join("chr1.net");
    assert!(chr1_path.exists());
    let chr1_content = std::fs::read_to_string(chr1_path)?;
    assert!(chr1_content.contains("net chr1 1000"));
    assert!(chr1_content.contains("fill 0 100 chrA + 0 100 id 1 score 100 ali 100"));

    let chr2_path = out_dir.path().join("chr2.net");
    assert!(chr2_path.exists());
    let chr2_content = std::fs::read_to_string(chr2_path)?;
    assert!(chr2_content.contains("net chr2 2000"));
    assert!(chr2_content.contains("fill 100 200 chrB - 100 200 id 2 score 200 ali 200"));

    Ok(())
}

// --- net subset tests ---

#[test]
fn test_net_subset_basic() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = cargo_bin_cmd!("pgr");

    // Chain: 0-1000.
    let mut chain_file = NamedTempFile::new()?;
    writeln!(
        chain_file,
        "chain 1000 chr1 1000 + 0 1000 chr2 1000 + 0 1000 1"
    )?;
    writeln!(chain_file, "1000")?;
    writeln!(chain_file)?;

    // Net: Fill 100-200.
    let mut net_file = NamedTempFile::new()?;
    writeln!(net_file, "net chr1 1000")?;
    writeln!(
        net_file,
        " fill 100 100 chr2 + 100 100 id 1 score 100 ali 100"
    )?;

    let out_file = NamedTempFile::new()?;
    let out_path = out_file.path().to_str().unwrap();

    cmd.arg("net")
        .arg("subset")
        .arg(net_file.path().to_str().unwrap())
        .arg(chain_file.path().to_str().unwrap())
        .arg(out_path)
        .assert()
        .success();

    let output = std::fs::read_to_string(out_path)?;
    println!("Output:\n{}", output);

    // Should contain a chain covering 100-200.
    // Header: size 1000. tStart 100. tEnd 200. qStart 100. qEnd 200.
    assert!(output.contains("chain 1000 chr1 1000 + 100 200 chr2 1000 + 100 200 1"));
    assert!(output.contains("100")); // Block size

    Ok(())
}

#[test]
fn test_net_subset_split_on_insert() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = cargo_bin_cmd!("pgr");

    // Chain: 0-1000.
    let mut chain_file = NamedTempFile::new()?;
    writeln!(
        chain_file,
        "chain 1000 chr1 1000 + 0 1000 chr2 1000 + 0 1000 1"
    )?;
    writeln!(chain_file, "1000")?;
    writeln!(chain_file)?;

    // Net: Fill 0-500.
    // Gap: 200-100 (at 200, len 100).
    //  Fill inside gap: 200-100 (id 2, but we only care about id 1 splitting)
    let mut net_file = NamedTempFile::new()?;
    writeln!(net_file, "net chr1 1000")?;
    writeln!(net_file, " fill 0 500 chr2 + 0 500 id 1 score 500 ali 500")?;
    writeln!(net_file, "  gap 200 100 chr2 + 200 100")?;
    writeln!(
        net_file,
        "   fill 200 100 chr2 + 200 100 id 2 score 100 ali 100"
    )?;

    let out_file = NamedTempFile::new()?;
    let out_path = out_file.path().to_str().unwrap();

    cmd.arg("net")
        .arg("subset")
        .arg(net_file.path().to_str().unwrap())
        .arg(chain_file.path().to_str().unwrap())
        .arg(out_path)
        .arg("--split-on-insert")
        .assert()
        .success();

    let output = std::fs::read_to_string(out_path)?;
    println!("Output:\n{}", output);

    // Should contain two parts of chain 1:
    // 1. 0-200.
    // 2. 300-500. (Gap was 200-300).

    assert!(output.contains("chain 1000 chr1 1000 + 0 200"));
    assert!(output.contains("chain 1000 chr1 1000 + 300 500"));

    Ok(())
}

// --- net filter tests ---

#[test]
fn test_net_filter_basic() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = cargo_bin_cmd!("pgr");

    // Create input net file
    let mut in_file = NamedTempFile::new()?;
    writeln!(in_file, "net chr1 1000")?;
    writeln!(in_file, " fill 0 100 chr2 + 0 100 id 1 score 100 ali 100")?; // Pass
    writeln!(
        in_file,
        " fill 200 100 chr2 + 200 100 id 2 score 50 ali 100"
    )?; // Fail score

    cmd.arg("net")
        .arg("filter")
        .arg(in_file.path().to_str().unwrap())
        .arg("--min-score")
        .arg("80")
        .assert()
        .success()
        .stdout(predicates::str::contains("id 1"))
        .stdout(predicates::str::contains("id 2").not());

    Ok(())
}

#[test]
fn test_net_filter_nested() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = cargo_bin_cmd!("pgr");

    // Nested structure
    let mut in_file = NamedTempFile::new()?;
    writeln!(in_file, "net chr1 1000")?;
    writeln!(in_file, " fill 0 500 chr2 + 0 500 id 1 score 200 ali 200")?; // Pass
    writeln!(in_file, "  gap 100 100 chr2 + 100 100")?;
    writeln!(
        in_file,
        "   fill 100 100 chr2 + 100 100 id 2 score 50 ali 50"
    )?; // Fail score

    // If child fails, parent should still be kept (pruning removes children).

    cmd.arg("net")
        .arg("filter")
        .arg(in_file.path().to_str().unwrap())
        .arg("--min-score")
        .arg("100")
        .assert()
        .success()
        .stdout(predicates::str::contains("id 1"))
        .stdout(predicates::str::contains("id 2").not());

    Ok(())
}

// --- net class tests ---

#[test]
fn test_net_class_basic() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = cargo_bin_cmd!("pgr");

    // Create input net file
    let mut in_file = NamedTempFile::new()?;
    writeln!(in_file, "net chr1 1000")?;
    writeln!(in_file, " fill 0 100 chr2 + 0 100 id 1 score 100 ali 100")?; // Unknown class
    writeln!(
        in_file,
        " fill 200 100 chr2 + 200 100 id 2 score 50 ali 100 type top"
    )?; // Class "top"
    writeln!(
        in_file,
        " fill 400 100 chr2 + 400 100 id 3 score 50 ali 100 type syn"
    )?; // Class "syn"

    // Total size 1000.
    // Fills: 100 (unknown) + 100 (top) + 100 (syn) = 300 bases.
    // Gaps: 1000 - 300 = 700 bases.

    cmd.arg("net")
        .arg("class")
        .arg(in_file.path().to_str().unwrap())
        .assert()
        .success()
        .stdout(predicates::str::contains("unknown").and(predicates::str::contains("100")))
        .stdout(predicates::str::contains("top").and(predicates::str::contains("100")))
        .stdout(predicates::str::contains("syn").and(predicates::str::contains("100")))
        .stdout(predicates::str::contains("gap").and(predicates::str::contains("700")));

    Ok(())
}
