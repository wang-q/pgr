use assert_cmd::Command;
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

fn get_path(filename: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("tests/pgr");
    path.push(filename);
    path
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
                    parts[12] = "ID";   // Ignore ID
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
fn test_2bit_size_pseudocat() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("pseudocat.2bit");
    let expected_output = get_path("pseudocat.sizes");
    let output = temp.path().join("out.sizes");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("size")
        .arg(&input)
        .arg("-o")
        .arg(&output);
    
    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    let expected_content = fs::read_to_string(&expected_output)?;
    
    assert_eq!(output_content.replace("\r\n", "\n"), expected_content.replace("\r\n", "\n"));

    Ok(())
}

#[test]
fn test_2bit_size_pseudopig() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("pseudopig.2bit");
    let expected_output = get_path("pseudopig.sizes");
    let output = temp.path().join("out.sizes");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("2bit")
        .arg("size")
        .arg(&input)
        .arg("-o")
        .arg(&output);
    
    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    let expected_content = fs::read_to_string(&expected_output)?;
    
    assert_eq!(output_content.replace("\r\n", "\n"), expected_content.replace("\r\n", "\n"));

    Ok(())
}

#[test]
fn test_lav_to_psl_lastz() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("lastz.lav");
    let expected_output = get_path("lastz.psl");
    let output = temp.path().join("out.psl");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("lav")
        .arg("to-psl")
        .arg(&input)
        .arg("-o")
        .arg(&output);
    
    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    let expected_content = fs::read_to_string(&expected_output)?;

    let output_lines: Vec<String> = output_content
        .lines()
        .filter(|l| !l.starts_with("#") && !l.trim().is_empty())
        .map(|l| l.to_string())
        .collect();

    let expected_lines: Vec<String> = expected_content
        .lines()
        .filter(|l| !l.starts_with("#") && !l.trim().is_empty())
        .map(|l| l.to_string())
        .collect();

    assert_eq!(output_lines.len(), expected_lines.len(), "Line count mismatch");

    for (i, (out, exp)) in output_lines.iter().zip(expected_lines.iter()).enumerate() {
        assert_eq!(out, exp, "Mismatch at line {}", i + 1);
    }

    Ok(())
}

#[test]
fn test_chaining_psl_lastz() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("lastz.psl");
    let t_2bit = get_path("pseudocat.2bit");
    let q_2bit = get_path("pseudopig.2bit");
    let expected_output = get_path("pslChain/lastz.raw.chain");
    let score_scheme = get_path("lastz.mat");
    let output = temp.path().join("out.chain");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("chaining")
        .arg("psl")
        .arg(&t_2bit)
        .arg(&q_2bit)
        .arg(&input)
        .arg("--output")
        .arg(&output)
        .arg("--score-scheme")
        .arg(&score_scheme)
        .arg("--linear-gap")
        .arg("medium")
        .arg("--min-score")
        .arg("1000"); // Matches axtChain -minScore=1000

    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    let expected_content = fs::read_to_string(&expected_output)?;

    let output_norm = normalize_chain_output(&output_content);
    let expected_norm = normalize_chain_output(&expected_content);

    // Verify chain lines are similar (ignoring score and ID)
    // And block data is present
    
    // Note: The order of chains might differ if not sorted. 
    // But pgr chaining psl usually outputs in input order or sorted by score/pos.
    // Let's check if we can match normalized content directly or if we need more robust comparison.
    // For now, try direct comparison of normalized strings.
    
    // If direct comparison fails due to ordering, we might need to sort chains.
    // But typically chaining follows input PSL order or specific logic.
    
    // Given the complexity of exact floating point scores or ID generation,
    // we primarily check if the structure and coordinates match.
    
    assert_eq!(output_norm.lines().count(), expected_norm.lines().count(), "Line count mismatch");
    
    // assert_eq!(output_norm, expected_norm); 
    // Commented out exact match for now to see if it runs. 
    // We might need to handle chain IDs which are usually sequential integers.
    // I added ID masking to normalize function.

    assert_eq!(output_norm, expected_norm);

    Ok(())
}
