use assert_cmd::Command;
use std::fs;
use std::io::Write;
use tempfile::{tempdir, NamedTempFile};

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
fn test_chain_net_reverse() -> Result<(), Box<dyn std::error::Error>> {
    let dir = tempdir()?;
    let chain_path = dir.path().join("in.chain");
    let t_sizes_path = dir.path().join("t.sizes");
    let q_sizes_path = dir.path().join("q.sizes");
    let t_net_path = dir.path().join("t.net");
    let q_net_path = dir.path().join("q.net");

    // Chain with reverse strand on query
    // chain 1000 chr1 1000 + 100 200 chr2 1000 - 800 900 1
    // 100
    // Q Coords on - strand: 800-900.
    // Mapped to + strand: (1000-900)=100, (1000-800)=200.
    // So Q net should show fill at 100, size 100.
    let chain_content = "chain 1000 chr1 1000 + 100 200 chr2 1000 - 800 900 1\n100\n\n";
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
    // T fill: 100 100. oChrom=chr2, oStrand=-, oStart=800 -> reverse -> 100, oSize=100.
    assert!(t_net_content.contains("fill 100 100 chr2 - 100 100"));

    let q_net_content = fs::read_to_string(&q_net_path)?;
    println!("Q Net content:\n{}", q_net_content);
    assert!(q_net_content.contains("net chr2 1000"));
    assert!(q_net_content.contains("fill 100 100 chr1 + 100 100"));
    
    Ok(())
}

// --- chain pre-net tests ---

#[test]
fn test_chain_pre_net() -> Result<(), Box<dyn std::error::Error>> {
    let temp = tempdir()?;
    let input = temp.path().join("input.chain");
    let t_sizes = temp.path().join("target.sizes");
    let q_sizes = temp.path().join("query.sizes");
    let output = temp.path().join("output.chain");

    fs::write(&t_sizes, "chr1 1000\n")?;
    fs::write(&q_sizes, "q1 1000\n")?;

    // Chain 1: 0-100 (Score 1000)
    // Chain 2: 0-100 (Score 500) -> Covered by 1, should be dropped
    // Chain 3: 50-150 (Score 400) -> Partially new (100-150), should be kept
    // Chain 4: 200-300 (Score 300) -> New, kept
    let chain_content = "\
chain 1000 chr1 1000 + 0 100 q1 1000 + 0 100 1
100
chain 500 chr1 1000 + 0 100 q1 1000 + 0 100 2
100
chain 400 chr1 1000 + 50 150 q1 1000 + 50 150 3
100
chain 300 chr1 1000 + 200 300 q1 1000 + 200 300 4
100
";
    fs::write(&input, chain_content)?;

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("chain")
        .arg("pre-net")
        .arg(&input)
        .arg(&t_sizes)
        .arg(&q_sizes)
        .arg(&output);

    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    
    // Check IDs
    // Should contain 1, 3, 4
    // Should NOT contain 2
    
    assert!(output_content.contains("chain 1000")); // ID 1
    assert!(!output_content.contains("chain 500"));  // ID 2
    assert!(output_content.contains("chain 400"));  // ID 3
    assert!(output_content.contains("chain 300"));  // ID 4

    Ok(())
}

#[test]
fn test_chain_pre_net_haplotype() -> Result<(), Box<dyn std::error::Error>> {
    let temp = tempdir()?;
    let input = temp.path().join("input.chain");
    let t_sizes = temp.path().join("target.sizes");
    let q_sizes = temp.path().join("query.sizes");
    let output = temp.path().join("output.chain");

    fs::write(&t_sizes, "chr1 1000\n")?;
    fs::write(&q_sizes, "q1_hap 1000\n")?;

    // Chain 1: q1_hap, score 1000
    let chain_content = "\
chain 1000 chr1 1000 + 0 100 q1_hap 1000 + 0 100 1
100
";
    fs::write(&input, chain_content)?;

    // Default: exclude hap
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("chain")
        .arg("pre-net")
        .arg(&input)
        .arg(&t_sizes)
        .arg(&q_sizes)
        .arg(&output);

    cmd.assert().success();
    let output_content = fs::read_to_string(&output)?;
    assert!(output_content.is_empty());

    // With --incl-hap
    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("chain")
        .arg("pre-net")
        .arg(&input)
        .arg(&t_sizes)
        .arg(&q_sizes)
        .arg(&output)
        .arg("--incl-hap");

    cmd.assert().success();
    let output_content = fs::read_to_string(&output)?;
    assert!(output_content.contains("chain 1000"));

    Ok(())
}

// --- chain sort tests ---

#[test]
fn test_chain_sort_basic() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("pgr")?;
    
    // Create unsorted chain file
    // Chain 1: score 100
    // Chain 2: score 200
    let mut input = NamedTempFile::new()?;
    writeln!(input, "chain 100 chr1 100 + 0 10 chr1 100 + 0 10 1\n10\n")?;
    writeln!(input, "chain 200 chr1 100 + 10 20 chr1 100 + 10 20 2\n10\n")?;
    
    let output = NamedTempFile::new()?;
    let out_path = output.path().to_str().unwrap();

    cmd.arg("chain")
        .arg("sort")
        .arg(input.path().to_str().unwrap())
        .arg("-o")
        .arg(out_path)
        .assert()
        .success();

    let content = std::fs::read_to_string(out_path)?;
    let lines: Vec<&str> = content.lines().collect();
    
    // Should see score 200 first, then 100
    // And IDs should be renumbered to 1, 2 (but original were 1, 2, so let's make sure order is correct)
    // First chain line should be score 200, ID 1 (renumbered)
    // Second chain line should be score 100, ID 2 (renumbered)
    
    // Find lines starting with "chain"
    let chains: Vec<&str> = lines.iter().filter(|l| l.starts_with("chain")).cloned().collect();
    assert_eq!(chains.len(), 2);
    
    assert!(chains[0].contains("200"));
    assert!(chains[0].ends_with("1")); // ID 1
    
    assert!(chains[1].contains("100"));
    assert!(chains[1].ends_with("2")); // ID 2
    
    Ok(())
}

#[test]
fn test_chain_sort_multiple_files() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("pgr")?;
    
    let mut f1 = NamedTempFile::new()?;
    writeln!(f1, "chain 100 chr1 100 + 0 10 chr1 100 + 0 10 1\n10\n")?;
    
    let mut f2 = NamedTempFile::new()?;
    writeln!(f2, "chain 300 chr1 100 + 20 30 chr1 100 + 20 30 2\n10\n")?;
    
    let output = NamedTempFile::new()?;
    let out_path = output.path().to_str().unwrap();

    cmd.arg("chain")
        .arg("sort")
        .arg(f1.path().to_str().unwrap())
        .arg(f2.path().to_str().unwrap())
        .arg("-o")
        .arg(out_path)
        .assert()
        .success();

    let content = std::fs::read_to_string(out_path)?;
    let chains: Vec<&str> = content.lines().filter(|l| l.starts_with("chain")).collect();
    
    // Order: 300, 100
    assert_eq!(chains.len(), 2);
    assert!(chains[0].contains("300"));
    assert!(chains[1].contains("100"));
    
    Ok(())
}

#[test]
fn test_chain_sort_save_id() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("pgr")?;
    
    let mut input = NamedTempFile::new()?;
    // Chain 1: score 100, ID 50
    // Chain 2: score 200, ID 60
    writeln!(input, "chain 100 chr1 100 + 0 10 chr1 100 + 0 10 50\n10\n")?;
    writeln!(input, "chain 200 chr1 100 + 10 20 chr1 100 + 10 20 60\n10\n")?;
    
    let output = NamedTempFile::new()?;
    let out_path = output.path().to_str().unwrap();

    cmd.arg("chain")
        .arg("sort")
        .arg(input.path().to_str().unwrap())
        .arg("-o")
        .arg(out_path)
        .arg("--save-id")
        .assert()
        .success();

    let content = std::fs::read_to_string(out_path)?;
    let chains: Vec<&str> = content.lines().filter(|l| l.starts_with("chain")).collect();
    
    // Order: 200, 100
    // IDs: 60, 50
    assert!(chains[0].contains("200"));
    assert!(chains[0].ends_with("60"));
    
    assert!(chains[1].contains("100"));
    assert!(chains[1].ends_with("50"));
    
    Ok(())
}

// --- chain split tests ---

#[test]
fn test_chain_split() -> Result<(), Box<dyn std::error::Error>> {
    let temp = tempdir()?;
    let input = temp.path().join("input.chain");
    let out_dir = temp.path().join("out");

    let chain_content = "\
chain 1000 chr1 1000 + 0 1000 q1 1000 + 0 1000 1
1000
chain 2000 chr2 1000 + 0 1000 q2 1000 + 0 1000 2
1000
";
    fs::write(&input, chain_content)?;

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("chain")
        .arg("split")
        .arg(&out_dir)
        .arg(&input);

    cmd.assert().success();

    assert!(out_dir.join("chr1.chain").exists());
    assert!(out_dir.join("chr2.chain").exists());

    let c1 = fs::read_to_string(out_dir.join("chr1.chain"))?;
    assert!(c1.contains("chr1"));
    assert!(!c1.contains("chr2"));

    let c2 = fs::read_to_string(out_dir.join("chr2.chain"))?;
    assert!(c2.contains("chr2"));
    assert!(!c2.contains("chr1"));

    Ok(())
}

#[test]
fn test_chain_split_lump() -> Result<(), Box<dyn std::error::Error>> {
    let temp = tempdir()?;
    let input = temp.path().join("input.chain");
    let out_dir = temp.path().join("out_lump");

    let chain_content = "\
chain 1000 chr1 1000 + 0 1000 q1 1000 + 0 1000 1
1000
chain 2000 chr2 1000 + 0 1000 q2 1000 + 0 1000 2
1000
";
    fs::write(&input, chain_content)?;

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("chain")
        .arg("split")
        .arg(&out_dir)
        .arg(&input)
        .arg("--lump")
        .arg("1");

    cmd.assert().success();

    assert!(out_dir.join("000.chain").exists());
    
    let content = fs::read_to_string(out_dir.join("000.chain"))?;
    assert!(content.contains("chr1"));
    assert!(content.contains("chr2"));

    Ok(())
}

// --- chain stitch tests ---

#[test]
fn test_chain_stitch() -> Result<(), Box<dyn std::error::Error>> {
    let temp = tempdir()?;
    let input = temp.path().join("fragments.chain");
    let output = temp.path().join("stitched.chain");

    // Create fragments with same ID but different ranges
    // Fragment 1: 0-1000
    // Fragment 2: 2000-3000
    let chain_content = "\
chain 1000 chr1 10000 + 0 1000 q1 10000 + 0 1000 1
1000
chain 1000 chr1 10000 + 2000 3000 q1 10000 + 2000 3000 1
1000
";
    fs::write(&input, chain_content)?;

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("chain")
        .arg("stitch")
        .arg(&input)
        .arg(&output);

    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    
    // Check if stitched into one chain
    // Should have:
    // header range: 0-3000
    // blocks: 1000, gap 1000, 1000
    // score: 2000
    
    assert!(output_content.contains("chain 2000 chr1 10000 + 0 3000 q1 10000 + 0 3000 1"));
    assert!(output_content.contains("1000 1000 1000"));
    assert!(output_content.lines().filter(|l| l.starts_with("chain")).count() == 1);

    Ok(())
}
