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

// 1. Alignment - lavToPsl
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

    let normalize = |s: &str| -> String {
        s.lines()
            .filter(|line| !line.starts_with('#')) // specific to psl header
            .collect::<Vec<&str>>()
            .join("\n")
    };

    assert_eq!(normalize(&output_content), normalize(&expected_content));

    Ok(())
}

// 1. Alignment (Prep) - fa size / faToTwoBit
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

// 2. Chain - axtChain
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

    assert_eq!(
        output_norm.lines().count(),
        expected_norm.lines().count(),
        "Line count mismatch"
    );

    assert_eq!(output_norm, expected_norm);

    Ok(())
}

// 2. Chain - chainAntiRepeat
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

    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    let expected_content = fs::read_to_string(&expected_output)?;

    let output_norm = normalize_chain_output(&output_content);
    let expected_norm = normalize_chain_output(&expected_content);

    assert_eq!(output_norm, expected_norm);

    Ok(())
}

// 2. Chain - chainMergeSort
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

// 2. Chain - chainPreNet
#[test]
fn test_chain_pre_net_lastz() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("all.chain");
    let t_sizes = get_path("pseudocat.sizes");
    let q_sizes = get_path("pseudopig.sizes");
    let expected_output = get_path("all.pre.chain");
    let output = temp.path().join("out.chain");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("chain")
        .arg("pre-net")
        .arg(&input)
        .arg(&t_sizes)
        .arg(&q_sizes)
        .arg(&output);

    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    let expected_content = fs::read_to_string(&expected_output)?;

    let output_norm = normalize_chain_output(&output_content);
    let expected_norm = normalize_chain_output(&expected_content);

    assert_eq!(output_norm, expected_norm);

    Ok(())
}

// 2. Chain - chainNet
#[test]
fn test_chain_net_lastz() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("all.pre.chain");
    let t_sizes = get_path("pseudocat.sizes");
    let q_sizes = get_path("pseudopig.sizes");
    let expected_t_net = get_path("pseudocat.chainnet");
    let expected_q_net = get_path("pseudopig.chainnet");

    let output_t_net = temp.path().join("out.t.net");
    let output_q_net = temp.path().join("out.q.net");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("chain")
                .arg("net")
                .arg(&input)
                .arg(&t_sizes)
                .arg(&q_sizes)
                .arg(&output_t_net)
                .arg(&output_q_net)
                .arg("--min-space")
                .arg("1")
                .arg("--min-score")
                .arg("2000");

    let output = cmd.output()?;
    let stderr = String::from_utf8_lossy(&output.stderr);
    println!("STDERR: {}", stderr);
    
    if !output.status.success() {
        panic!("Command failed with status: {}", output.status);
    }

    let t_net_content = fs::read_to_string(&output_t_net)?;
    let expected_t_content = fs::read_to_string(&expected_t_net)?;
    assert_eq!(t_net_content, expected_t_content);

    let q_net_content = fs::read_to_string(&output_q_net)?;
    let expected_q_content = fs::read_to_string(&expected_q_net)?;
    assert_eq!(q_net_content, expected_q_content);

    Ok(())
}

// 2. Chain - netSyntenic
#[test]
fn test_net_syntenic_lastz() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("pseudocat.chainnet");
    let expected_output = get_path("noClass.net");
    let output = temp.path().join("out.net");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("net").arg("syntenic").arg(&input).arg(&output);

    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    let expected_content = fs::read_to_string(&expected_output)?;
    assert_eq!(output_content, expected_content);

    Ok(())
}

// 2. Chain - netChainSubset
#[test]
fn test_net_subset_lastz() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let net_input = get_path("noClass.net");
    let chain_input = get_path("all.chain");
    let expected_output = get_path("subset.chain");
    let output = temp.path().join("out.chain");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("net")
        .arg("subset")
        .arg(&net_input)
        .arg(&chain_input)
        .arg(&output);

    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    let expected_content = fs::read_to_string(&expected_output)?;

    let output_norm = normalize_chain_output(&output_content);
    let expected_norm = normalize_chain_output(&expected_content);

    assert_eq!(output_norm, expected_norm);

    Ok(())
}

// 2. Chain - chainStitchId
#[test]
fn test_chain_stitch_lastz() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("subset.chain");
    let expected_output = get_path("over.chain");
    let output = temp.path().join("out.chain");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("chain").arg("stitch").arg(&input).arg(&output);

    cmd.assert().success();

    let output_content = fs::read_to_string(&output)?;
    let expected_content = fs::read_to_string(&expected_output)?;

    let output_norm = normalize_chain_output(&output_content);
    let expected_norm = normalize_chain_output(&expected_content);

    assert_eq!(output_norm, expected_norm);

    Ok(())
}

fn normalize_net_output(content: &str) -> String {
    content
        .lines()
        .filter(|line| !line.starts_with('#'))
        .collect::<Vec<&str>>()
        .join("\n")
}

#[test]
fn test_net_split_lastz() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let input = get_path("noClass.net");
    let output_dir = temp.path().join("net");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("net").arg("split").arg(&input).arg(&output_dir);

    cmd.assert().success();

    let output_file = output_dir.join("cat.net");
    assert!(output_file.exists());

    let output_content = fs::read_to_string(&output_file)?;
    let expected_output = get_path("net/cat.net");
    let expected_content = fs::read_to_string(&expected_output)?;

    assert_eq!(
        normalize_net_output(&output_content),
        normalize_net_output(&expected_content)
    );

    Ok(())
}

// 3. Axt - netToAxt
#[test]
fn test_net_to_axt_lastz() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let net_input = get_path("net/cat.net");
    let chain_input = get_path("all.pre.chain");
    let t_2bit = get_path("pseudocat.2bit");
    let q_2bit = get_path("pseudopig.2bit");
    let output = temp.path().join("cat.axt");

    let mut cmd = Command::cargo_bin("pgr")?;
    cmd.arg("net")
        .arg("to-axt")
        .arg(&net_input)
        .arg(&chain_input)
        .arg(&t_2bit)
        .arg(&q_2bit)
        .arg(&output);

    cmd.assert().success();

    assert!(output.exists());
    assert!(fs::metadata(&output)?.len() > 0);

    let expected_output = get_path("axtNet/cat.axt");
    let output_content = fs::read_to_string(&output)?;
    let expected_content = fs::read_to_string(&expected_output)?;
    
    // Normalize newlines
    let output_norm = output_content.replace("\r\n", "\n");
    let expected_norm = expected_content.replace("\r\n", "\n");

    assert_eq!(output_norm, expected_norm);

    Ok(())
}

// 3. Axt - axtSort
#[test]
fn test_axt_sort_lastz() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    
    // Use the expected output (already sorted) as input to verify idempotency/format handling
    // This avoids dependency on net-to-axt command in this test
    let input_path = get_path("axtNet/cat.axt");
    let output = temp.path().join("cat.axt");

    let mut cmd_sort = Command::cargo_bin("pgr")?;
    cmd_sort.arg("axt")
        .arg("sort")
        .arg(&input_path)
        .arg("-o")
        .arg(&output);

    cmd_sort.assert().success();

    let output_content = fs::read_to_string(&output)?;
    let expected_content = fs::read_to_string(&input_path)?;
    
    // Normalize newlines for cross-platform comparison
    let output_norm = output_content.replace("\r\n", "\n");
    let expected_norm = expected_content.replace("\r\n", "\n");

    assert_eq!(output_norm, expected_norm);

    Ok(())
}
