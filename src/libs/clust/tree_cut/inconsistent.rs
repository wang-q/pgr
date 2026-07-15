use super::{assign_clusters, Partition};
use crate::libs::phylo::tree::Tree;
use anyhow::Result;
use std::collections::HashMap;

/// Cut tree based on inconsistent coefficient threshold.
pub fn cut_inconsistent(tree: &Tree, threshold: f64, depth: usize) -> Result<Partition> {
    let root = tree
        .get_root()
        .ok_or_else(|| anyhow::anyhow!("Tree has no root"))?;

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::phylo::tree::Tree;
    use std::collections::HashMap;

    fn parse_tree(nwk: &str) -> Tree {
        Tree::from_newick(nwk).expect("valid newick")
    }

    fn cluster_names(part: &Partition, tree: &Tree) -> Vec<Vec<String>> {
        let mut groups: HashMap<usize, Vec<String>> = HashMap::new();
        for (&leaf_id, &cid) in &part.assignment {
            let name = tree
                .get_node(leaf_id)
                .and_then(|n| n.name.clone())
                .unwrap_or_else(|| format!("Node_{}", leaf_id));
            groups.entry(cid).or_default().push(name);
        }
        let mut clusters: Vec<Vec<String>> = groups.into_values().collect();
        for c in &mut clusters {
            c.sort();
        }
        clusters.sort();
        clusters
    }

    #[test]
    fn test_cut_inconsistent_low_threshold_splits() {
        // Tree: ((A:1,B:1):1,C:1);
        // Root has positive inconsistency, so threshold=0 forces a split.
        let tree = parse_tree("((A:1,B:1):1,C:1);");
        let part = cut_inconsistent(&tree, 0.0, 2).unwrap();
        let mut clusters = cluster_names(&part, &tree);
        clusters.sort();
        assert_eq!(clusters, vec![vec!["A", "B"], vec!["C"]]);
    }

    #[test]
    fn test_cut_inconsistent_high_threshold_keeps_root() {
        let tree = parse_tree("((A:1,B:1):1,C:1);");
        let part = cut_inconsistent(&tree, 100.0, 2).unwrap();
        let clusters = cluster_names(&part, &tree);
        assert_eq!(clusters, vec![vec!["A", "B", "C"]]);
    }

    #[test]
    fn test_cut_inconsistent_depth_changes_result() {
        // Tree: (((A:1,B:1):1,C:1):1,D:1);
        // Depth 1 vs depth 2 gives different root inconsistency, so the same
        // threshold can produce different numbers of clusters.
        let tree = parse_tree("(((A:1,B:1):1,C:1):1,D:1);");

        let part_depth1 = cut_inconsistent(&tree, 0.8, 1).unwrap();
        let clusters_depth1 = cluster_names(&part_depth1, &tree);
        assert_eq!(clusters_depth1, vec![vec!["A", "B", "C", "D"]]);

        let part_depth2 = cut_inconsistent(&tree, 0.8, 2).unwrap();
        let clusters_depth2 = cluster_names(&part_depth2, &tree);
        assert_eq!(clusters_depth2, vec![vec!["A", "B", "C"], vec!["D"]]);
    }

    #[test]
    fn test_cut_inconsistent_empty_tree() {
        let tree = Tree::new();
        assert!(cut_inconsistent(&tree, 1.0, 2).is_err());
    }
}
