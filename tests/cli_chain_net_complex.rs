use assert_cmd::cargo::cargo_bin_cmd;
use std::fs;
use tempfile::tempdir;

#[test]
fn test_chain_net_greedy_overlap() -> Result<(), Box<dyn std::error::Error>> {
    let dir = tempdir()?;
    let chain_path = dir.path().join("in.chain");
    let t_sizes_path = dir.path().join("t.sizes");
    let q_sizes_path = dir.path().join("q.sizes");
    let t_net_path = dir.path().join("t.net");
    let q_net_path = dir.path().join("q.net");

    // Target: chr1 2000
    fs::write(&t_sizes_path, "chr1 2000\n")?;
    fs::write(&q_sizes_path, "chr2 2000\n")?;

    // Chain A: Score 2000. 0-1000.
    // chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
    let c1 = "chain 2000 chr1 2000 + 0 1000 chr2 2000 + 0 1000 1\n1000\n\n";

    // Chain B: Score 1000. 500-1500.
    // Overlaps 500-1000 with Chain A.
    let c2 = "chain 1000 chr1 2000 + 500 1500 chr2 2000 + 500 1500 2\n1000\n\n";

    // Write chains in order (A then B). Since A has higher score, order shouldn't matter if we sort,
    // but pgr net might expect sorted input or sort internally.
    // The current pgr implementation checks if sorted but doesn't enforce it strictly?
    // Let's provide them sorted by score (2000 then 1000).
    fs::write(&chain_path, format!("{}{}", c1, c2))?;

    let mut cmd = cargo_bin_cmd!("pgr");
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

    let t_net_content = fs::read_to_string(&t_net_path)?;
    println!("T Net content:\n{}", t_net_content);

    // Expectation:
    // Chain 1 (Score 2000) takes 0-1000.
    // Chain 2 (Score 1000) takes 1000-1500 (500-1000 is blocked).

    // Check for Chain 1 fill
    assert!(t_net_content.contains("fill 0 1000 chr2 + 0 1000 id 1"));

    // Check for Chain 2 fill
    // It should start at 1000, size 500.
    // Corresponding Q coords: Chain 2 maps T:500-1500 to Q:500-1500.
    // So T:1000-1500 maps to Q:1000-1500.
    assert!(t_net_content.contains("fill 1000 500 chr2 + 1000 500 id 2"));

    Ok(())
}

#[test]
fn test_chain_net_nested_fill() -> Result<(), Box<dyn std::error::Error>> {
    let dir = tempdir()?;
    let chain_path = dir.path().join("in.chain");
    let t_sizes_path = dir.path().join("t.sizes");
    let q_sizes_path = dir.path().join("q.sizes");
    let t_net_path = dir.path().join("t.net");
    let q_net_path = dir.path().join("q.net");

    fs::write(&t_sizes_path, "chr1 2000\n")?;
    fs::write(&q_sizes_path, "chr2 2000\n")?;

    // Chain A: Score 2000. 0-1000.
    // Has a gap:
    // Block 1: 0-400 (400bp)
    // Gap: 400-600 (200bp) on T.
    // Block 2: 600-1000 (400bp).
    // chain header...
    // 400 200 200 (size dt dq) -> dt=200, dq=200
    // 400
    let c1 = "chain 2000 chr1 2000 + 0 1000 chr2 2000 + 0 1000 1\n400 200 200\n400\n\n";

    // Chain B: Score 1000. 450-550.
    // Fits inside the gap of Chain A (400-600).
    // chain ... 450 550 ...
    let c2 = "chain 1000 chr1 2000 + 450 550 chr2 2000 + 450 550 2\n100\n\n";

    fs::write(&chain_path, format!("{}{}", c1, c2))?;

    let mut cmd = cargo_bin_cmd!("pgr");
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

    let t_net_content = fs::read_to_string(&t_net_path)?;
    println!("T Net content:\n{}", t_net_content);

    // Expectation:
    // Level 1: Chain 1 (Fill 0-1000).
    // Level 2: Inside Gap of Chain 1 (Gap 400-200), we see Chain 2 (Fill 450-100).

    // Check Chain 1
    assert!(t_net_content.contains("fill 0 1000 chr2 + 0 1000 id 1"));

    // Check Chain 2 (nested)
    // Should be indented (checked by structure, but here we just check presence first)
    // fill 450 100 chr2 + 450 100 id 2
    assert!(t_net_content.contains("fill 450 100 chr2 + 450 100 id 2"));

    Ok(())
}
