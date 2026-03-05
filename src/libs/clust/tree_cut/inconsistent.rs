use super::{assign_clusters, Partition};
use crate::libs::phylo::tree::Tree;
use std::collections::HashMap;

/// Cut tree based on inconsistent coefficient threshold.
pub fn cut_inconsistent(tree: &Tree, threshold: f64, depth: usize) -> Result<Partition, String> {
    let root = tree.get_root().ok_or("Tree has no root")?;

    // 1. Compute node heights
    let heights = crate::libs::phylo::tree::stat::compute_node_heights(tree);

    // 2. Calculate inconsistency
    let inconsistencies =
        crate::libs::phylo::tree::stat::calculate_inconsistency(tree, &heights, depth);

    // 3. Compute max inconsistency in subtree (post-order)
    let mut max_inc_subtree = HashMap::new();
    let post_order = crate::libs::phylo::tree::traversal::postorder(tree, root);

    for id in post_order {
        if let Some(node) = tree.get_node(id) {
            let my_inc = *inconsistencies.get(&id).unwrap_or(&0.0);
            let mut max_inc = my_inc;

            for &child in &node.children {
                let child_max = *max_inc_subtree.get(&child).unwrap_or(&0.0);
                if child_max > max_inc {
                    max_inc = child_max;
                }
            }
            max_inc_subtree.insert(id, max_inc);
        }
    }

    // 4. Top-down traversal to find clusters
    let mut clusters = Vec::new();
    let mut stack = vec![root];

    while let Some(node_id) = stack.pop() {
        let max_inc = *max_inc_subtree.get(&node_id).unwrap_or(&0.0);

        if max_inc <= threshold {
            // This node and all descendants satisfy the condition
            clusters.push(node_id);
        } else {
            // Violation somewhere in subtree.
            // Check children.
            if let Some(node) = tree.get_node(node_id) {
                if node.children.is_empty() {
                    // Leaf has max_inc=0 usually.
                    clusters.push(node_id);
                } else {
                    for &child in &node.children {
                        stack.push(child);
                    }
                }
            }
        }
    }

    assign_clusters(tree, clusters)
}
