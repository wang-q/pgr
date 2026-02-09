use assert_cmd::Command;
use pgr::libs::phylo::tree::Tree;

// Helper to check support values and branch lengths with delta
fn check_support_and_length(
    output: &str,
    expected_support: &str,
    expected_length: f64,
    epsilon: f64,
) -> bool {
    // Parse the output Newick string into a Tree object
    let trees = Tree::from_newick_multi(output).expect("Failed to parse output Newick");
    
    // Iterate through all nodes to find the matching pattern
    // We look for a node where:
    // 1. Name matches expected_support
    // 2. Length matches expected_length within epsilon
    for tree in &trees {
        // Use postorder traversal to iterate all nodes since nodes field is private
        if let Some(root) = tree.get_root() {
            if let Ok(nodes) = tree.postorder(&root) {
                for node_id in nodes {
                    if let Some(node) = tree.get_node(node_id) {
                        if let Some(name) = &node.name {
                            if name == expected_support {
                                if let Some(len) = node.length {
                                     let diff: f64 = (len - expected_length).abs();
                                     if diff < epsilon {
                                         return true;
                                     }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    false
}

#[test]
fn test_nwk_support_simple() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("nwk")
        .arg("support")
        .arg("tests/newick/HRV.nwk")
        .arg("tests/newick/HRV_20reps.nwk")
        .output()?;
        
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout)?;

    // Check for some known support values from test_nw_support_simple.exp
    // e.g., ...HRV1B_1:0.123339)6:0.076821... (support 6, length 0.076821)
    // ...HRV64_1:0.064173)16:0... (support 16, length 0.0)
    
    let epsilon = 1e-6;
    
    assert!(check_support_and_length(&stdout, "6", 0.076821, epsilon), "Failed to find node with support 6 and length 0.076821");
    assert!(check_support_and_length(&stdout, "16", 0.0, epsilon), "Failed to find node with support 16 and length 0.0");
    assert!(check_support_and_length(&stdout, "20", 0.227116, epsilon), "Failed to find node with support 20 and length 0.227116");

    Ok(())
}

#[test]
fn test_nwk_support_percent() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("nwk")
        .arg("support")
        .arg("tests/newick/HRV.nwk")
        .arg("tests/newick/HRV_20reps.nwk")
        .arg("--percent")
        .output()?;

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout)?;

    // Check for some known support values from test_nw_support_percent.exp
    // 6/20 * 100 = 30
    // 16/20 * 100 = 80
    // 20/20 * 100 = 100
    
    let epsilon = 1e-6;

    assert!(check_support_and_length(&stdout, "30", 0.076821, epsilon), "Failed to find node with support 30 and length 0.076821");
    assert!(check_support_and_length(&stdout, "80", 0.0, epsilon), "Failed to find node with support 80 and length 0.0");
    assert!(check_support_and_length(&stdout, "100", 0.227116, epsilon), "Failed to find node with support 100 and length 0.227116");

    Ok(())
}

#[test]
fn test_nwk_support_multi() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("pgr")?;
    let output = cmd.arg("nwk")
        .arg("support")
        .arg("tests/newick/3_HRV.nwk")
        .arg("tests/newick/HRV_20reps.nwk")
        .output()?;

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout)?;
    
    let trees = Tree::from_newick_multi(&stdout).expect("Failed to parse output Newick");
    assert_eq!(trees.len(), 3);
    
    let epsilon = 1e-6;

    // Check first tree
    assert!(check_support_and_length(&stdout, "6", 0.076821, epsilon), "Failed to find node in multi-tree output");

    Ok(())
}
