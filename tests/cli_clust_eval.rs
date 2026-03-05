use assert_cmd::Command;

#[test]
fn test_clust_eval_perfect_match() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("clust")
        .arg("eval")
        .arg("tests/clust/perfect_1.tsv")
        .arg("tests/clust/perfect_2.tsv")
        .arg("--format")
        .arg("cluster")
        .output()?;

    let stdout = String::from_utf8(output.stdout)?;
    assert!(output.status.success());

    // Check header
    assert!(stdout.contains("ari\tami\thomogeneity\tcompleteness\tv_measure"));

    // Check values (should be 1.0)
    // The second line should be numbers
    let lines: Vec<&str> = stdout.lines().collect();
    let values: Vec<&str> = lines[1].split_whitespace().collect();

    assert_eq!(values.len(), 5);
    // ARI
    assert!((values[0].parse::<f64>()? - 1.0).abs() < 1e-6);
    // AMI
    assert!((values[1].parse::<f64>()? - 1.0).abs() < 1e-6);
    // V-Measure
    assert!((values[4].parse::<f64>()? - 1.0).abs() < 1e-6);

    Ok(())
}

#[test]
fn test_clust_eval_disjoint() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("clust")
        .arg("eval")
        .arg("tests/clust/perfect_1.tsv")
        .arg("tests/clust/disjoint_2.tsv")
        .arg("--format")
        .arg("cluster")
        .output()?;

    let stdout = String::from_utf8(output.stdout)?;
    assert!(output.status.success());

    let lines: Vec<&str> = stdout.lines().collect();
    let values: Vec<&str> = lines[1].split_whitespace().collect();

    // ARI should be -0.5
    assert!((values[0].parse::<f64>()? - (-0.5)).abs() < 1e-6);
    // AMI should be -0.5 (approx, due to log calculation precision?)
    // In my manual calculation, it was exactly -0.5.
    // Let's check roughly.
    let ami = values[1].parse::<f64>()?;
    assert!((ami - (-0.5)).abs() < 1e-4, "AMI was {}", ami);

    Ok(())
}

#[test]
fn test_clust_eval_single_vs_singletons() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("clust")
        .arg("eval")
        .arg("tests/clust/single_1.tsv")
        .arg("tests/clust/singletons.tsv")
        .arg("--format")
        .arg("cluster")
        .output()?;

    let stdout = String::from_utf8(output.stdout)?;
    assert!(output.status.success());

    let lines: Vec<&str> = stdout.lines().collect();
    let values: Vec<&str> = lines[1].split_whitespace().collect();

    // ARI = 0
    assert!((values[0].parse::<f64>()? - 0.0).abs() < 1e-6);
    // AMI = 0
    assert!((values[1].parse::<f64>()? - 0.0).abs() < 1e-6);

    Ok(())
}

#[test]
fn test_clust_eval_pair_format() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("clust")
        .arg("eval")
        .arg("tests/clust/pair_1.tsv")
        .arg("tests/clust/pair_2.tsv")
        .output()?;

    let stdout = String::from_utf8(output.stdout)?;
    assert!(output.status.success());

    // Check header
    assert!(stdout.contains("ari\tami\thomogeneity\tcompleteness\tv_measure"));

    // Check values (should be 1.0)
    let lines: Vec<&str> = stdout.lines().collect();
    let values: Vec<&str> = lines[1].split_whitespace().collect();

    assert_eq!(values.len(), 5);
    // ARI
    assert!((values[0].parse::<f64>()? - 1.0).abs() < 1e-6);
    // AMI
    assert!((values[1].parse::<f64>()? - 1.0).abs() < 1e-6);

    Ok(())
}

#[test]
fn test_clust_eval_internal_silhouette() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("clust")
        .arg("eval")
        .arg("tests/clust/eval/simple.pair")
        .arg("--matrix")
        .arg("tests/clust/eval/simple.matrix.phy")
        .output()?;

    let stdout = String::from_utf8(output.stdout)?;
    assert!(output.status.success());

    // Output format:
    // silhouette
    // 0.5167

    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines[0], "silhouette");

    let score = lines[1].parse::<f64>()?;
    let expected = 1.55 / 3.0;
    assert!(
        (score - expected).abs() < 1e-4,
        "Score was {}, expected {}",
        score,
        expected
    );

    Ok(())
}
