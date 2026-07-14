#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;
use tempfile::tempdir;

// --- chain net tests ---

#[test]
fn test_chain_net_basic() {
    let dir = tempdir().unwrap();
    let chain_path = dir.path().join("in.chain");
    let t_sizes_path = dir.path().join("t.sizes");
    let q_sizes_path = dir.path().join("q.sizes");
    let t_net_path = dir.path().join("t.net");
    let q_net_path = dir.path().join("q.net");

    // Create inputs
    // Chain with one block
    // chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
    let chain_content = "chain 1000 chr1 1000 + 0 100 chr2 1000 + 0 100 1\n100\n\n";
    fs::write(&chain_path, chain_content).unwrap();

    fs::write(&t_sizes_path, "chr1 1000\n").unwrap();
    fs::write(&q_sizes_path, "chr2 1000\n").unwrap();

    PgrCmd::new()
        .args(&[
            "chain",
            "net",
            chain_path.to_str().unwrap(),
            t_sizes_path.to_str().unwrap(),
            q_sizes_path.to_str().unwrap(),
            t_net_path.to_str().unwrap(),
            q_net_path.to_str().unwrap(),
            "--min-score=0",
            "--min-space=1",
        ])
        .run();

    // Verify output
    let t_net_content = fs::read_to_string(&t_net_path).unwrap();
    println!("T Net content:\n{}", t_net_content);
    assert!(t_net_content.contains("net chr1 1000"));
    // fill start size oChrom oStrand oStart oSize id score ali
    assert!(t_net_content.contains("fill 0 100 chr2 + 0 100"));

    let q_net_content = fs::read_to_string(&q_net_path).unwrap();
    println!("Q Net content:\n{}", q_net_content);
    assert!(q_net_content.contains("net chr2 1000"));
    assert!(q_net_content.contains("fill 0 100 chr1 + 0 100"));
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

    let mut cmd_2bit = assert_cmd::Command::cargo_bin("pgr").unwrap();
    cmd_2bit
        .arg("fa")
        .arg("to-2bit")
        .arg(t_fa_path.to_str().unwrap())
        .arg("-o")
        .arg(t_2bit_path.to_str().unwrap());
    cmd_2bit.assert().success();

    let mut cmd_2bit = assert_cmd::Command::cargo_bin("pgr").unwrap();
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

    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    cmd.arg("chain")
        .arg("anti-repeat")
        .arg("--target-2bit")
        .arg(t_2bit_path.to_str().unwrap())
        .arg("--query-2bit")
        .arg(q_2bit_path.to_str().unwrap())
        .arg(chain_path.to_str().unwrap())
        .arg("-o")
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
fn test_chain_sort() {
    let dir = tempdir().unwrap();
    let chain1_path = dir.path().join("1.chain");
    let chain2_path = dir.path().join("2.chain");
    let out_path = dir.path().join("out.chain");

    // Chain 1: score 100
    // chain 100 chr1 100 + 0 10 chr2 100 + 0 10 1
    let c1 = "chain 100 chr1 100 + 0 10 chr2 100 + 0 10 1\n10\n\n";
    fs::write(&chain1_path, c1).unwrap();

    // Chain 2: score 200
    // chain 200 chr1 100 + 20 30 chr2 100 + 20 30 2
    let c2 = "chain 200 chr1 100 + 20 30 chr2 100 + 20 30 2\n10\n\n";
    fs::write(&chain2_path, c2).unwrap();

    PgrCmd::new()
        .args(&[
            "chain",
            "sort",
            chain1_path.to_str().unwrap(),
            chain2_path.to_str().unwrap(),
            "--outfile",
            out_path.to_str().unwrap(),
        ])
        .run();

    let output = fs::read_to_string(&out_path).unwrap();
    // Should be sorted by score descending: 200 then 100
    // IDs are renumbered by default: 1 then 2

    let lines: Vec<&str> = output.lines().filter(|l| l.starts_with("chain")).collect();
    assert_eq!(lines.len(), 2);
    assert!(lines[0].contains("chain 200"));
    assert!(lines[1].contains("chain 100"));
}

#[test]
fn test_chain_sort_input_list() {
    let dir = tempdir().unwrap();
    let chain1_path = dir.path().join("1.chain");
    let chain2_path = dir.path().join("2.chain");
    let list_path = dir.path().join("list.txt");
    let out_path = dir.path().join("out.chain");

    // Chain 1: score 100
    let c1 = "chain 100 chr1 100 + 0 10 chr2 100 + 0 10 1\n10\n\n";
    fs::write(&chain1_path, c1).unwrap();

    // Chain 2: score 200
    let c2 = "chain 200 chr1 100 + 20 30 chr2 100 + 20 30 2\n10\n\n";
    fs::write(&chain2_path, c2).unwrap();

    // Create list file
    // Note: On Windows, paths might contain backslashes.
    // The implementation uses simple line reading, which should handle standard paths fine.
    let list_content = format!(
        "{}\n{}",
        chain1_path.to_str().unwrap(),
        chain2_path.to_str().unwrap()
    );
    fs::write(&list_path, list_content).unwrap();

    PgrCmd::new()
        .args(&[
            "chain",
            "sort",
            "--input-list",
            list_path.to_str().unwrap(),
            "--outfile",
            out_path.to_str().unwrap(),
        ])
        .run();

    let output = fs::read_to_string(&out_path).unwrap();
    let lines: Vec<&str> = output.lines().filter(|l| l.starts_with("chain")).collect();
    assert_eq!(lines.len(), 2);
    assert!(lines[0].contains("chain 200"));
    assert!(lines[1].contains("chain 100"));
}

#[test]
fn test_chain_sort_mixed_inputs() {
    let dir = tempdir().unwrap();
    let chain1_path = dir.path().join("1.chain");
    let chain2_path = dir.path().join("2.chain");
    let list_path = dir.path().join("list.txt");
    let out_path = dir.path().join("out.chain");

    // Chain 1: score 100
    let c1 = "chain 100 chr1 100 + 0 10 chr2 100 + 0 10 1\n10\n\n";
    fs::write(&chain1_path, c1).unwrap();

    // Chain 2: score 200
    let c2 = "chain 200 chr1 100 + 20 30 chr2 100 + 20 30 2\n10\n\n";
    fs::write(&chain2_path, c2).unwrap();

    // Create list file with only chain1
    let list_content = format!("{}\n", chain1_path.to_str().unwrap());
    fs::write(&list_path, list_content).unwrap();

    PgrCmd::new()
        .args(&[
            "chain",
            "sort",
            "--input-list",
            list_path.to_str().unwrap(),
            chain2_path.to_str().unwrap(), // Pass chain2 as arg
            "--outfile",
            out_path.to_str().unwrap(),
        ])
        .run();

    let output = fs::read_to_string(&out_path).unwrap();
    let lines: Vec<&str> = output.lines().filter(|l| l.starts_with("chain")).collect();
    assert_eq!(lines.len(), 2);
    assert!(lines[0].contains("chain 200"));
    assert!(lines[1].contains("chain 100"));
}

// --- chain split tests ---

#[test]
fn test_chain_split_by_target() {
    let dir = tempdir().unwrap();
    let in1_path = dir.path().join("in1.chain");
    let in2_path = dir.path().join("in2.chain");
    let outdir = dir.path().join("split_out");

    // in1.chain: chr1 (score 100) + chr2 (score 200)
    let c1 = "chain 100 chr1 1000 + 0 10 chr2 1000 + 0 10 1\n10\n\n";
    let c2 = "chain 200 chr2 1000 + 20 30 chr3 1000 + 20 30 2\n10\n\n";
    fs::write(&in1_path, format!("{}{}", c1, c2)).unwrap();

    // in2.chain: chr1 (score 300)
    let c3 = "chain 300 chr1 1000 + 40 50 chr2 1000 + 40 50 3\n10\n\n";
    fs::write(&in2_path, c3).unwrap();

    PgrCmd::new()
        .args(&[
            "chain",
            "split",
            in1_path.to_str().unwrap(),
            in2_path.to_str().unwrap(),
            "--outdir",
            outdir.to_str().unwrap(),
        ])
        .run();

    // Verify chr1.chain contains 2 chains (from in1 and in2)
    let chr1_path = outdir.join("chr1.chain");
    assert!(chr1_path.exists(), "chr1.chain should exist");
    let chr1_content = fs::read_to_string(&chr1_path).unwrap();
    let chr1_chains: Vec<&str> = chr1_content
        .lines()
        .filter(|l| l.starts_with("chain"))
        .collect();
    assert_eq!(chr1_chains.len(), 2);
    assert!(chr1_chains[0].contains("chain 100"));
    assert!(chr1_chains[1].contains("chain 300"));

    // Verify chr2.chain contains 1 chain
    let chr2_path = outdir.join("chr2.chain");
    assert!(chr2_path.exists(), "chr2.chain should exist");
    let chr2_content = fs::read_to_string(&chr2_path).unwrap();
    let chr2_chains: Vec<&str> = chr2_content
        .lines()
        .filter(|l| l.starts_with("chain"))
        .collect();
    assert_eq!(chr2_chains.len(), 1);
    assert!(chr2_chains[0].contains("chain 200"));
}

#[test]
fn test_chain_split_by_query() {
    let dir = tempdir().unwrap();
    let in1_path = dir.path().join("in1.chain");
    let outdir = dir.path().join("split_out");

    // Two chains with different query names
    let c1 = "chain 100 chr1 1000 + 0 10 chr_qA 1000 + 0 10 1\n10\n\n";
    let c2 = "chain 200 chr1 1000 + 20 30 chr_qB 1000 + 20 30 2\n10\n\n";
    fs::write(&in1_path, format!("{}{}", c1, c2)).unwrap();

    PgrCmd::new()
        .args(&[
            "chain",
            "split",
            in1_path.to_str().unwrap(),
            "--by-query",
            "--outdir",
            outdir.to_str().unwrap(),
        ])
        .run();

    // Verify chr_qA.chain contains 1 chain
    let qa_path = outdir.join("chr_qA.chain");
    assert!(qa_path.exists(), "chr_qA.chain should exist");
    let qa_content = fs::read_to_string(&qa_path).unwrap();
    let qa_chains: Vec<&str> = qa_content
        .lines()
        .filter(|l| l.starts_with("chain"))
        .collect();
    assert_eq!(qa_chains.len(), 1);
    assert!(qa_chains[0].contains("chain 100"));

    // Verify chr_qB.chain contains 1 chain
    let qb_path = outdir.join("chr_qB.chain");
    assert!(qb_path.exists(), "chr_qB.chain should exist");
    let qb_content = fs::read_to_string(&qb_path).unwrap();
    let qb_chains: Vec<&str> = qb_content
        .lines()
        .filter(|l| l.starts_with("chain"))
        .collect();
    assert_eq!(qb_chains.len(), 1);
    assert!(qb_chains[0].contains("chain 200"));
}

// --- chain net / pre-net sort-order tests ---

#[test]
fn test_chain_net_unsorted_fails() {
    let dir = tempdir().unwrap();
    let chain_path = dir.path().join("in.chain");
    let t_sizes_path = dir.path().join("t.sizes");
    let q_sizes_path = dir.path().join("q.sizes");
    let t_net_path = dir.path().join("t.net");
    let q_net_path = dir.path().join("q.net");

    // Scores are ascending, not descending.
    let c1 = "chain 100 chr1 1000 + 0 10 chr2 1000 + 0 10 1\n10\n\n";
    let c2 = "chain 200 chr1 1000 + 20 30 chr2 1000 + 20 30 2\n10\n\n";
    fs::write(&chain_path, format!("{}{}", c1, c2)).unwrap();

    fs::write(&t_sizes_path, "chr1 1000\n").unwrap();
    fs::write(&q_sizes_path, "chr2 1000\n").unwrap();

    let (_stdout, stderr) = PgrCmd::new()
        .args(&[
            "chain",
            "net",
            chain_path.to_str().unwrap(),
            t_sizes_path.to_str().unwrap(),
            q_sizes_path.to_str().unwrap(),
            t_net_path.to_str().unwrap(),
            q_net_path.to_str().unwrap(),
        ])
        .run_fail();

    assert!(
        stderr.contains("Input not sorted by score"),
        "expected sort error in stderr, got: {}",
        stderr
    );
}

#[test]
fn test_chain_pre_net_unsorted_fails() {
    let dir = tempdir().unwrap();
    let chain_path = dir.path().join("in.chain");
    let t_sizes_path = dir.path().join("t.sizes");
    let q_sizes_path = dir.path().join("q.sizes");
    let out_path = dir.path().join("out.chain");

    // Scores are ascending, not descending.
    let c1 = "chain 100 chr1 1000 + 0 10 chr2 1000 + 0 10 1\n10\n\n";
    let c2 = "chain 200 chr1 1000 + 20 30 chr2 1000 + 20 30 2\n10\n\n";
    fs::write(&chain_path, format!("{}{}", c1, c2)).unwrap();

    fs::write(&t_sizes_path, "chr1 1000\n").unwrap();
    fs::write(&q_sizes_path, "chr2 1000\n").unwrap();

    let (_stdout, stderr) = PgrCmd::new()
        .args(&[
            "chain",
            "pre-net",
            chain_path.to_str().unwrap(),
            t_sizes_path.to_str().unwrap(),
            q_sizes_path.to_str().unwrap(),
            "-o",
            out_path.to_str().unwrap(),
        ])
        .run_fail();

    assert!(
        stderr.contains("Input not sorted by score"),
        "expected sort error in stderr, got: {}",
        stderr
    );
}

// --- chain stitch tests ---

#[test]
fn test_chain_stitch() {
    let dir = tempdir().unwrap();
    let chain_path = dir.path().join("in.chain");
    let out_path = dir.path().join("out.chain");

    // Two fragments with the same ID; scores should be summed.
    let c1 = "chain 100 chr1 1000 + 0 10 chr2 1000 + 0 10 1\n10\n\n";
    let c2 = "chain 200 chr1 1000 + 20 30 chr2 1000 + 20 30 1\n10\n\n";
    fs::write(&chain_path, format!("{}{}", c1, c2)).unwrap();

    PgrCmd::new()
        .args(&[
            "chain",
            "stitch",
            chain_path.to_str().unwrap(),
            "-o",
            out_path.to_str().unwrap(),
        ])
        .run();

    let output = fs::read_to_string(&out_path).unwrap();
    let lines: Vec<&str> = output.lines().filter(|l| l.starts_with("chain")).collect();
    assert_eq!(lines.len(), 1);
    assert!(lines[0].contains("chain 300"));
}

// --- chain split lump tests ---

#[test]
fn test_chain_split_lump() {
    let dir = tempdir().unwrap();
    let in1_path = dir.path().join("in1.chain");
    let outdir = dir.path().join("split_out");

    // chr10 -> 10 % 10 = 0, chr21 -> 21 % 10 = 1.
    let c1 = "chain 100 chr10 1000 + 0 10 chr2 1000 + 0 10 1\n10\n\n";
    let c2 = "chain 200 chr21 1000 + 20 30 chr2 1000 + 20 30 2\n10\n\n";
    fs::write(&in1_path, format!("{}{}", c1, c2)).unwrap();

    PgrCmd::new()
        .args(&[
            "chain",
            "split",
            in1_path.to_str().unwrap(),
            "--outdir",
            outdir.to_str().unwrap(),
            "--lump",
            "10",
        ])
        .run();

    let p0 = outdir.join("000.chain");
    let p1 = outdir.join("001.chain");
    assert!(p0.exists(), "000.chain should exist");
    assert!(p1.exists(), "001.chain should exist");

    let content0 = fs::read_to_string(&p0).unwrap();
    assert!(content0.contains("chain 100"));

    let content1 = fs::read_to_string(&p1).unwrap();
    assert!(content1.contains("chain 200"));
}

// --- chain anti-repeat negative strand test ---

#[test]
fn test_chain_anti_repeat_negative_strand() -> Result<(), Box<dyn std::error::Error>> {
    let dir = tempdir()?;
    let chain_path = dir.path().join("in.chain");
    let out_path = dir.path().join("out.chain");
    let t_fa_path = dir.path().join("target.fa");
    let q_fa_path = dir.path().join("query.fa");

    let t_seq = ">chr1\nACTGACTGAC\n";
    fs::write(&t_fa_path, t_seq)?;
    // Query is the reverse complement of target so the negative-strand chain is valid.
    let q_seq = ">chr2\nGTCAGTCAGT\n";
    fs::write(&q_fa_path, q_seq)?;

    let t_2bit_path = dir.path().join("target.2bit");
    let q_2bit_path = dir.path().join("query.2bit");

    let mut cmd_2bit = assert_cmd::Command::cargo_bin("pgr").unwrap();
    cmd_2bit
        .arg("fa")
        .arg("to-2bit")
        .arg(t_fa_path.to_str().unwrap())
        .arg("-o")
        .arg(t_2bit_path.to_str().unwrap());
    cmd_2bit.assert().success();

    let mut cmd_2bit = assert_cmd::Command::cargo_bin("pgr").unwrap();
    cmd_2bit
        .arg("fa")
        .arg("to-2bit")
        .arg(q_fa_path.to_str().unwrap())
        .arg("-o")
        .arg(q_2bit_path.to_str().unwrap());
    cmd_2bit.assert().success();

    // Query on the negative strand.
    let chain = "chain 1000 chr1 10 + 0 10 chr2 10 - 0 10 1\n10\n\n";
    fs::write(&chain_path, chain)?;

    let mut cmd = assert_cmd::Command::cargo_bin("pgr").unwrap();
    cmd.arg("chain")
        .arg("anti-repeat")
        .arg("--target-2bit")
        .arg(t_2bit_path.to_str().unwrap())
        .arg("--query-2bit")
        .arg(q_2bit_path.to_str().unwrap())
        .arg(chain_path.to_str().unwrap())
        .arg("-o")
        .arg(out_path.to_str().unwrap())
        .arg("--min-score")
        .arg("800");

    cmd.assert().success();

    let output = fs::read_to_string(&out_path)?;
    let lines: Vec<&str> = output.lines().filter(|l| l.starts_with("chain")).collect();
    assert_eq!(lines.len(), 1);
    assert!(lines[0].contains("chain 1000"));
    assert!(lines[0].contains("chr2 10 - 0 10 1"));

    Ok(())
}
