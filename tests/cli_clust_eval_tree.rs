use assert_cmd::Command;
use std::fs;
use tempfile::tempdir;

#[test]
fn test_clust_eval_tree_internal() -> anyhow::Result<()> {
    let temp = tempdir()?;
    let partition_path = temp.path().join("partition.tsv");
    let tree_path = "tests/newick/dist.nwk";

    // Cluster 1: A, B
    // Cluster 2: D, E
    // C and F are in tree but not in partition (should be ignored)
    // Format: ClusterID <tab> Item
    let partition_content = "1\tA\n1\tB\n2\tD\n2\tE\n";
    fs::write(&partition_path, partition_content)?;

    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd
        .arg("clust")
        .arg("eval")
        .arg(partition_path.to_str().unwrap())
        .arg("--tree")
        .arg(tree_path)
        .output()?;

    if !output.status.success() {
        eprintln!("STDERR:\n{}", String::from_utf8_lossy(&output.stderr));
    }
    assert!(output.status.success());

    let stdout = String::from_utf8(output.stdout)?;
    let lines: Vec<&str> = stdout.lines().collect();
    
    // Check Header
    assert_eq!(lines[0], "silhouette\tdunn\tc_index\tgamma\ttau");
    
    // Values
    let values: Vec<&str> = lines[1].split('\t').collect();
    
    // Silhouette: 0.494949
    let s_score = values[0].parse::<f64>()?;
    let expected_s = 49.0 / 99.0;
    assert!((s_score - expected_s).abs() < 1e-4, "Silhouette was {}, expected {}", s_score, expected_s);

    // Dunn: 1.3333
    let d_score = values[1].parse::<f64>()?;
    let expected_d = 4.0 / 3.0;
    assert!((d_score - expected_d).abs() < 1e-4, "Dunn was {}, expected {}", d_score, expected_d);

    // C-Index: 0.0
    let c_score = values[2].parse::<f64>()?;
    assert!((c_score - 0.0).abs() < 1e-4, "C-index was {}", c_score);

    // Gamma: ~0.877
    let g_score = values[3].parse::<f64>()?;
    assert!((g_score - 0.877058).abs() < 1e-4, "Gamma was {}", g_score);

    // Tau: ~0.756
    let t_score = values[4].parse::<f64>()?;
    assert!((t_score - 0.755928).abs() < 1e-4, "Tau was {}", t_score);

    Ok(())
}
