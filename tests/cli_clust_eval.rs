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
    if !output.status.success() {
        eprintln!("STDERR:\n{}", String::from_utf8_lossy(&output.stderr));
    }
    assert!(output.status.success());

    // Check header
    assert!(stdout.contains("ari\tami\thomogeneity\tcompleteness\tv_measure\tfmi\tnmi\tmi"));

    // Check values (should be 1.0)
    // The second line should be numbers
    let lines: Vec<&str> = stdout.lines().collect();
    let values: Vec<&str> = lines[1].split_whitespace().collect();

    assert_eq!(values.len(), 8);
    // ARI
    assert!((values[0].parse::<f64>()? - 1.0).abs() < 1e-6);
    // AMI
    assert!((values[1].parse::<f64>()? - 1.0).abs() < 1e-6);
    // V-Measure
    assert!((values[4].parse::<f64>()? - 1.0).abs() < 1e-6);
    // FMI
    assert!((values[5].parse::<f64>()? - 1.0).abs() < 1e-6);
    // NMI
    assert!((values[6].parse::<f64>()? - 1.0).abs() < 1e-6);

    Ok(())
}

#[test]
fn test_clust_eval_no_singletons() -> anyhow::Result<()> {
    // Create temporary files
    // Truth: {A, B, C} (Cluster 1), {D} (Singleton), {E} (Singleton)
    // Pred: {A, B, C} (Cluster 1), {D, E} (Cluster 2)
    // Format: ClusterID <tab> Item

    let truth_content = "1\tA\n1\tB\n1\tC\n2\tD\n3\tE\n";
    let pred_content = "1\tA\n1\tB\n1\tC\n2\tD\n2\tE\n";

    let truth_path = "tests/clust/eval_ns_truth.tsv";
    let pred_path = "tests/clust/eval_ns_pred.tsv";

    std::fs::write(truth_path, truth_content)?;
    std::fs::write(pred_path, pred_content)?;

    // 1. Run WITHOUT --no-singletons
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("clust")
        .arg("eval")
        .arg(pred_path)
        .arg(truth_path)
        .output()?;

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout)?;
    let lines: Vec<&str> = stdout.lines().collect();
    let values: Vec<&str> = lines[1].split_whitespace().collect();
    let ari = values[0].parse::<f64>()?;
    // ARI should be < 1.0
    assert!(ari < 0.99, "ARI was {} (expected < 0.99)", ari);

    // 2. Run WITH --no-singletons
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("clust")
        .arg("eval")
        .arg(pred_path)
        .arg(truth_path)
        .arg("--no-singletons")
        .output()?;

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout)?;
    let lines: Vec<&str> = stdout.lines().collect();
    let values: Vec<&str> = lines[1].split_whitespace().collect();
    let ari = values[0].parse::<f64>()?;

    // ARI should be 1.0
    assert!((ari - 1.0).abs() < 1e-6, "ARI was {}", ari);

    // Cleanup
    std::fs::remove_file(truth_path)?;
    std::fs::remove_file(pred_path)?;

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
    if !output.status.success() {
        eprintln!("STDERR:\n{}", String::from_utf8_lossy(&output.stderr));
    }
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
    if !output.status.success() {
        eprintln!("STDERR:\n{}", String::from_utf8_lossy(&output.stderr));
    }
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
    if !output.status.success() {
        eprintln!("STDERR:\n{}", String::from_utf8_lossy(&output.stderr));
    }
    assert!(output.status.success());

    // Check header
    assert!(stdout.contains("ari\tami\thomogeneity\tcompleteness\tv_measure\tfmi\tnmi\tmi"));

    // Check values (should be 1.0)
    let lines: Vec<&str> = stdout.lines().collect();
    let values: Vec<&str> = lines[1].split_whitespace().collect();

    assert_eq!(values.len(), 8);
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

#[test]
fn test_clust_eval_internal_db() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("clust")
        .arg("eval")
        .arg("tests/clust/eval/db.pair")
        .arg("--coords")
        .arg("tests/clust/eval/db.coords.tsv")
        .output()?;

    let stdout = String::from_utf8(output.stdout)?;
    assert!(output.status.success());

    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines[0], "davies_bouldin");

    let score = lines[1].parse::<f64>()?;
    let expected = 0.2;
    assert!(
        (score - expected).abs() < 1e-4,
        "Score was {}, expected {}",
        score,
        expected
    );

    Ok(())
}
