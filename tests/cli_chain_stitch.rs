use assert_cmd::Command;
use std::fs;
use tempfile::tempdir;

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
