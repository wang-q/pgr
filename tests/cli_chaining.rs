use assert_cmd::Command;
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

fn get_path(subcommand: &str, dir: &str, filename: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("tests/chaining");
    path.push(subcommand);
    path.push(dir);
    path.push(filename);
    path
}

fn create_2bit(dir: &TempDir, name: &str, content: &str) -> anyhow::Result<std::path::PathBuf> {
    let fa_path = dir.path().join(format!("{}.fa", name));
    let bit_path = dir.path().join(format!("{}.2bit", name));
    fs::write(&fa_path, content)?;

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("fa")
        .arg("to-2bit")
        .arg(&fa_path)
        .arg("-o")
        .arg(&bit_path);
    cmd.assert().success();

    Ok(bit_path)
}

#[test]
fn test_chaining_psl_basic() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    
    // Create genomes
    // Target: chr1: 1000bp
    // Query: chr2: 1000bp
    let t_seq = ">chr1\n".to_string() + &"A".repeat(1000);
    let q_seq = ">chr2\n".to_string() + &"A".repeat(1000);
    
    let t_2bit = create_2bit(&temp, "t", &t_seq)?;
    let q_2bit = create_2bit(&temp, "q", &q_seq)?;
    
    // Create PSL
    // match: 100, t: 0-100, q: 0-100
    // Use tabs for PSL fields as required by Psl::from_str
    // qName=chr2 (matches q.2bit), tName=chr1 (matches t.2bit)
    let psl_content = "100\t0\t0\t0\t0\t0\t0\t0\t+\tchr2\t1000\t0\t100\tchr1\t1000\t0\t100\t1\t100,\t0,\t0,\n";
    let psl_path = temp.path().join("in.psl");
    fs::write(&psl_path, psl_content)?;
    
    let output_path = temp.path().join("out.chain");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("chaining")
        .arg("psl")
        .arg(&t_2bit)
        .arg(&q_2bit)
        .arg(&psl_path)
        .arg("-o")
        .arg(&output_path)
        .arg("--min-score=0");

    let assert = cmd.assert().success();
    let stderr = String::from_utf8(assert.get_output().stderr.clone())?;
    println!("Stderr: {}", stderr);

    let output = fs::read_to_string(&output_path)?;
    println!("Chain output:\n{}", output);
    
    // chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
    assert!(output.contains("chain"));
    assert!(output.contains("chr1 1000 + 0 100 chr2 1000 + 0 100"));
    assert!(output.contains("100")); // block size

    Ok(())
}

fn normalize_chain_output(content: &str) -> String {
    content
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            if line.starts_with("chain") {
                let mut parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() > 2 {
                    parts[1] = "SCORE"; // Ignore score
                }
                parts.join(" ")
            } else {
                line.to_string()
            }
        })
        .collect::<Vec<String>>()
        .join("\n")
}

#[test]
fn test_chaining_psl_new_style_lastz() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    // Reuse data from tests/psl/chain
    let input = get_path("psl", "input", "newStyleLastz.psl");
    let t_2bit = get_path("psl", "input", "hg19.chrM.2bit");
    let q_2bit = get_path("psl", "input", "susScr3.chrM.2bit");
    let score_scheme = get_path("psl", "input", "newStyleLastz.Q.txt");
    let expected_output = get_path("psl", "expected", "newStyleLastz.chain");
    let output = temp.path().join("out.chain");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("chaining")
        .arg("psl")
        .arg(&t_2bit)
        .arg(&q_2bit)
        .arg(&input)
        .arg("--score-scheme")
        .arg(&score_scheme)
        .arg("--output")
        .arg(&output)
        .arg("--linear-gap")
        .arg("loose")
        .arg("--min-score")
        .arg("3000");
    
    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    let expected_content = fs::read_to_string(&expected_output)?;
    
    let output_norm = normalize_chain_output(&output_content);
    let expected_norm = normalize_chain_output(&expected_content);

    assert_eq!(output_norm, expected_norm);

    Ok(())
}
