use assert_cmd::Command;
use std::fs;
use tempfile::tempdir;

// --- chain net tests ---

#[test]
fn test_chain_net_basic() -> Result<(), Box<dyn std::error::Error>> {
    let dir = tempdir()?;
    let chain_path = dir.path().join("in.chain");
    let t_sizes_path = dir.path().join("t.sizes");
    let q_sizes_path = dir.path().join("q.sizes");
    let t_net_path = dir.path().join("t.net");
    let q_net_path = dir.path().join("q.net");

    // Create inputs
    // Chain with one block
    // chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
    let chain_content = "chain 1000 chr1 1000 + 0 100 chr2 1000 + 0 100 1\n100\n\n";
    fs::write(&chain_path, chain_content)?;

    fs::write(&t_sizes_path, "chr1 1000\n")?;
    fs::write(&q_sizes_path, "chr2 1000\n")?;

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("chain")
        .arg("net")
        .arg(chain_path.to_str().unwrap())
        .arg(t_sizes_path.to_str().unwrap())
        .arg(q_sizes_path.to_str().unwrap())
        .arg(t_net_path.to_str().unwrap())
        .arg(q_net_path.to_str().unwrap())
        .arg("--min-score=0")
        .arg("--min-space=1");

    cmd.assert().success();

    // Verify output
    let t_net_content = fs::read_to_string(&t_net_path)?;
    println!("T Net content:\n{}", t_net_content);
    assert!(t_net_content.contains("net chr1 1000"));
    // fill start size oChrom oStrand oStart oSize id score ali
    assert!(t_net_content.contains("fill 0 100 chr2 + 0 100"));

    let q_net_content = fs::read_to_string(&q_net_path)?;
    println!("Q Net content:\n{}", q_net_content);
    assert!(q_net_content.contains("net chr2 1000"));
    assert!(q_net_content.contains("fill 0 100 chr1 + 0 100"));

    Ok(())
}

#[test]
fn test_chain_anti_repeat() -> Result<(), Box<dyn std::error::Error>> {
    let dir = tempdir()?;
    let chain_path = dir.path().join("in.chain");
    let out_path = dir.path().join("out.chain");
    let t_fa_path = dir.path().join("target.fa");
    let q_fa_path = dir.path().join("query.fa");

    // Target: chr1
    // 0-10: "AAAAAAAAAA" (low complexity)
    // 10-20: "actgactgac" (repeats)
    // 20-30: "ACTGACTGAC" (good)
    let t_seq = ">chr1\nAAAAAAAAAAactgactgacACTGACTGAC\n";
    fs::write(&t_fa_path, t_seq)?;

    // Query: chr2 (matches perfectly)
    let q_seq = ">chr2\nAAAAAAAAAAactgactgacACTGACTGAC\n";
    fs::write(&q_fa_path, q_seq)?;

    // Convert FASTA to 2bit
    let t_2bit_path = dir.path().join("target.2bit");
    let q_2bit_path = dir.path().join("query.2bit");

    let mut cmd_2bit = Command::cargo_bin("pgr")?;
    cmd_2bit
        .arg("fa")
        .arg("to-2bit")
        .arg(t_fa_path.to_str().unwrap())
        .arg("-o")
        .arg(t_2bit_path.to_str().unwrap());
    cmd_2bit.assert().success();

    let mut cmd_2bit = Command::cargo_bin("pgr")?;
    cmd_2bit
        .arg("fa")
        .arg("to-2bit")
        .arg(q_fa_path.to_str().unwrap())
        .arg("-o")
        .arg(q_2bit_path.to_str().unwrap());
    cmd_2bit.assert().success();

    // Chain 1: 0-10 (Low complexity) -> score 1000
    // Chain 2: 10-20 (Repeat) -> score 1000
    // Chain 3: 20-30 (Good) -> score 1000

    // chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
    // Chain 1
    let c1 = "chain 1000 chr1 30 + 0 10 chr2 30 + 0 10 1\n10\n\n";
    // Chain 2
    let c2 = "chain 1000 chr1 30 + 10 20 chr2 30 + 10 20 2\n10\n\n";
    // Chain 3
    let c3 = "chain 1000 chr1 30 + 20 30 chr2 30 + 20 30 3\n10\n\n";

    fs::write(&chain_path, format!("{}{}{}", c1, c2, c3))?;

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("chain")
        .arg("anti-repeat")
        .arg("--target")
        .arg(t_2bit_path.to_str().unwrap())
        .arg("--query")
        .arg(q_2bit_path.to_str().unwrap())
        .arg(chain_path.to_str().unwrap())
        .arg(out_path.to_str().unwrap())
        .arg("--min-score")
        .arg("800");

    cmd.assert().success();

    let output = fs::read_to_string(&out_path)?;
    // Chain 1 should be filtered (score ~10)
    assert!(!output.contains(" 10 1\n"));
    // Chain 2 should be filtered (score 0)
    assert!(!output.contains(" 20 2\n"));
    // Chain 3 should be kept (score 2000)
    assert!(output.contains(" 30 3\n"));

    Ok(())
}

#[test]
fn test_chain_sort() -> Result<(), Box<dyn std::error::Error>> {
    let dir = tempdir()?;
    let chain1_path = dir.path().join("1.chain");
    let chain2_path = dir.path().join("2.chain");
    let out_path = dir.path().join("out.chain");

    // Chain 1: score 100
    // chain 100 chr1 100 + 0 10 chr2 100 + 0 10 1
    let c1 = "chain 100 chr1 100 + 0 10 chr2 100 + 0 10 1\n10\n\n";
    fs::write(&chain1_path, c1)?;

    // Chain 2: score 200
    // chain 200 chr1 100 + 20 30 chr2 100 + 20 30 2
    let c2 = "chain 200 chr1 100 + 20 30 chr2 100 + 20 30 2\n10\n\n";
    fs::write(&chain2_path, c2)?;

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("chain")
        .arg("sort")
        .arg(chain1_path.to_str().unwrap())
        .arg(chain2_path.to_str().unwrap())
        .arg("--output")
        .arg(out_path.to_str().unwrap());

    cmd.assert().success();

    let output = fs::read_to_string(&out_path)?;
    // Should be sorted by score descending: 200 then 100
    // IDs are renumbered by default: 1 then 2

    let lines: Vec<&str> = output.lines().filter(|l| l.starts_with("chain")).collect();
    assert_eq!(lines.len(), 2);
    assert!(lines[0].contains("chain 200"));
    assert!(lines[1].contains("chain 100"));

    Ok(())
}
