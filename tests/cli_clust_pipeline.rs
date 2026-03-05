mod common;
use crate::common::*;
use std::fs;

// --- Helper: Generate Synthetic Blobs ---
// Generates 3 well-separated clusters in 2D space.
// Format: ID \t X,Y
// Ground Truth: ID \t ClusterID
fn generate_blobs(data_file: &str, truth_file: &str) -> Result<(), Box<dyn std::error::Error>> {
    let mut data_out = String::new();
    let mut truth_out = String::new();

    // Cluster 1: Center (0,0), 10 points
    for i in 0..10 {
        let x = 0.0 + (i as f64 * 0.1);
        let y = 0.0 + (i as f64 * 0.1);
        let id = format!("C1_{}", i);
        data_out.push_str(&format!("{}\t{:.4},{:.4}\n", id, x, y));
        truth_out.push_str(&format!("1\t{}\n", id));
    }

    // Cluster 2: Center (10,10), 10 points
    for i in 0..10 {
        let x = 10.0 + (i as f64 * 0.1);
        let y = 10.0 + (i as f64 * 0.1);
        let id = format!("C2_{}", i);
        data_out.push_str(&format!("{}\t{:.4},{:.4}\n", id, x, y));
        truth_out.push_str(&format!("2\t{}\n", id));
    }

    // Cluster 3: Center (20,0), 10 points
    for i in 0..10 {
        let x = 20.0 + (i as f64 * 0.1);
        let y = 0.0 + (i as f64 * 0.1);
        let id = format!("C3_{}", i);
        data_out.push_str(&format!("{}\t{:.4},{:.4}\n", id, x, y));
        truth_out.push_str(&format!("3\t{}\n", id));
    }

    fs::create_dir_all("tests/pipeline")?;
    fs::write(data_file, data_out)?;
    fs::write(truth_file, truth_out)?;
    Ok(())
}

#[test]
fn test_clust_pipeline_full() {
    let temp_dir = tempfile::Builder::new()
        .prefix("pgr_pipeline_test")
        .tempdir()
        .expect("Failed to create temp dir");
    let base_dir = temp_dir.path().to_str().unwrap();

    let data_file = format!("{}/blobs.tsv", base_dir);
    let truth_file = format!("{}/blobs.truth.tsv", base_dir);

    // Intermediate files
    let dist_file = format!("{}/blobs.dist.tsv", base_dir);
    let phy_file = format!("{}/blobs.phy", base_dir);
    let tree_file = format!("{}/blobs.nwk", base_dir);
    let cut_file = format!("{}/blobs.cut.tsv", base_dir);
    let _eval_file = format!("{}/blobs.eval.tsv", base_dir);

    // 1. Generate Data
    generate_blobs(&data_file, &truth_file).expect("Failed to generate data");

    // 2. Calculate Distances (pgr dist vector)
    // Input: ID \t X,Y
    // Output: ID1 \t ID2 \t Dist
    let (_stdout, stderr) = PgrCmd::new()
        .args(&[
            "dist", "vector", &data_file, "--mode", "euclid", "-o", &dist_file,
        ])
        .run();
    assert!(stderr.is_empty(), "dist vector failed: {}", stderr);
    assert!(fs::metadata(&dist_file).is_ok(), "dist file not created");

    // 3. Convert to PHYLIP Matrix (pgr mat to-phylip)
    let (_stdout, stderr) = PgrCmd::new()
        .args(&["mat", "to-phylip", &dist_file, "-o", &phy_file])
        .run();
    assert!(stderr.is_empty(), "mat to-phylip failed: {}", stderr);
    assert!(fs::metadata(&phy_file).is_ok(), "phylip file not created");

    // 4. Hierarchical Clustering (pgr clust hier)
    // Method: Ward (standard for Euclidean)
    let (_stdout, stderr) = PgrCmd::new()
        .args(&[
            "clust", "hier", &phy_file, "--method", "ward", "-o",
            &tree_file, // Explicit output file
        ])
        .run();
    assert!(stderr.is_empty(), "clust hier failed: {}", stderr);
    assert!(fs::metadata(&tree_file).is_ok(), "tree file not created");

    // 5. Cut Tree (pgr nwk cut)
    // We know there are 3 clusters, so use --k 3
    let (_stdout, stderr) = PgrCmd::new()
        .args(&[
            "nwk", "cut", &tree_file, "--k", "3", "--format",
            "pair", // Output: Rep \t Member (compatible with eval)
            "-o", &cut_file,
        ])
        .run();
    assert!(stderr.is_empty(), "nwk cut failed: {}", stderr);
    assert!(fs::metadata(&cut_file).is_ok(), "cut file not created");

    // 6. Evaluate (pgr clust eval)
    // Compare cut result with ground truth
    // Output ARI should be 1.0 for perfect clustering
    let (stdout, stderr) = PgrCmd::new()
        .args(&[
            "clust",
            "eval",
            &truth_file, // Ground Truth
            &cut_file,   // Prediction
            "--format",
            "pair", // Both are in pair format (or at least compatible)
        ])
        .run();

    assert!(stderr.is_empty(), "clust eval failed: {}", stderr);

    // Evaluate output format:
    // Metric   Value
    // ari      1.0000
    // ...
    println!("Eval Output:\n{}", stdout);

    // Debug: Print cut file content if evaluation fails
    if !stdout.contains("1.000000") {
        // Check for perfect score
        let cut_content = fs::read_to_string(&cut_file).unwrap_or_default();
        println!("Cut File Content:\n{}", cut_content);
    }

    let ari_line = stdout
        .lines()
        .find(|l| l.trim().starts_with("0.") || l.trim().starts_with("1."))
        .expect("Score line not found");
    let parts: Vec<&str> = ari_line.split_whitespace().collect();
    let ari_val: f64 = parts[0].parse().expect("Failed to parse ARI"); // ARI is first column

    assert!(
        (ari_val - 1.0).abs() < 1e-4,
        "ARI should be 1.0, got {}",
        ari_val
    );

    // Cleanup handled by tempdir drop, or we can explicit close if needed but Drop trait handles it.
}
