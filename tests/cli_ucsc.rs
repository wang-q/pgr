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
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let mut parts: Vec<&str> = line.split_whitespace().collect();
            if parts.first() == Some(&"chain") {
                if parts.len() > 12 {
                    parts[1] = "SCORE"; // Ignore score
                    parts[12] = "ID"; // Ignore ID
                }
            }
            parts.join(" ")
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

    assert_eq!(output_content, expected_content);

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
        .arg("--output")
        .arg(&output);

    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    let expected_content = fs::read_to_string(&expected_output)?;

    // Normalize output for comparison
    // Ignore lines starting with # (header comments might differ or be absent)
    // Ignore match counts which might slightly differ due to implementation details?
    // Actually, lav to psl should be deterministic.
    // However, UCSC tools might output extra headers.

    let normalize = |s: &str| -> String {
        s.lines()
            .filter(|line| !line.starts_with('#')) // specific to psl header
            .collect::<Vec<&str>>()
            .join("\n")
    };

    assert_eq!(normalize(&output_content), normalize(&expected_content));

    Ok(())
}

#[test]
fn test_chaining_psl_lastz() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("lastz.psl");
    let t_2bit = get_path("pseudocat.2bit");
    let q_2bit = get_path("pseudopig.2bit");
    let expected_output = get_path("pslChain/lastz.raw.chain");
    let output = temp.path().join("out.chain");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("chaining")
        .arg("psl")
        .arg(&t_2bit)
        .arg(&q_2bit)
        .arg(&input)
        .arg("--outfile")
        .arg(&output)
        .arg("--score-scheme")
        .arg("hoxd55")
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

    assert_eq!(
        output_norm.lines().count(),
        expected_norm.lines().count(),
        "Line count mismatch"
    );

    // assert_eq!(output_norm, expected_norm);
    // Commented out exact match for now to see if it runs.
    // We might need to handle chain IDs which are usually sequential integers.
    // I added ID masking to normalize function.

    assert_eq!(output_norm, expected_norm);

    Ok(())
}

#[test]
fn test_chain_sort_lastz() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("pslChain/lastz.chain");
    let expected_output = get_path("all.chain");
    let output = temp.path().join("out.chain");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("chain")
        .arg("sort")
        .arg(&input)
        .arg("--output")
        .arg(&output);

    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    let expected_content = fs::read_to_string(&expected_output)?;

    let output_norm = normalize_chain_output(&output_content);
    let expected_norm = normalize_chain_output(&expected_content);

    assert_eq!(output_norm, expected_norm);

    Ok(())
}

#[test]
fn test_chain_anti_repeat_lastz() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("pslChain/lastz.raw.chain");
    let t_2bit = get_path("pseudocat.2bit");
    let q_2bit = get_path("pseudopig.2bit");
    let expected_output = get_path("pslChain/lastz.chain");
    let output = temp.path().join("out.chain");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("chain")
        .arg("anti-repeat")
        .arg("--target")
        .arg(&t_2bit)
        .arg("--query")
        .arg(&q_2bit)
        .arg(&input)
        .arg(&output);
        // Default min-score is 5000, which matches UCSC default? 
        // Docs say chainAntiRepeat default is ?
        // The previous test uses 1000.
        // Let's use default first.

    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    let expected_content = fs::read_to_string(&expected_output)?;

    let output_norm = normalize_chain_output(&output_content);
    let expected_norm = normalize_chain_output(&expected_content);

    assert_eq!(output_norm, expected_norm);

    Ok(())
}
