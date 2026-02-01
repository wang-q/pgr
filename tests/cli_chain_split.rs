use assert_cmd::Command;
use std::fs;
use tempfile::tempdir;

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
