use super::{assign_clusters, compute_heights, Partition};
use crate::libs::phylo::node::NodeId;
use crate::libs::phylo::tree::Tree;
use anyhow::Result;

/// Cut tree into K clusters
pub fn cut_k(tree: &Tree, k: usize) -> Result<Partition> {
    if k == 0 {
        anyhow::bail!("K must be >= 1");
    }

    let root = tree
        .get_root()
        .ok_or_else(|| anyhow::anyhow!("Tree has no root"))?;
    let leaves = crate::libs::phylo::tree::stat::get_leaves(tree, root);
    let n_leaves = leaves.len();

    if k >= n_leaves {
        // Return all leaves as individual clusters
        let mut part = Partition::new();
        for (i, leaf) in leaves.into_iter().enumerate() {
            part.assignment.insert(leaf, i + 1);
        }
        part.num_clusters = n_leaves;
        return Ok(part);
    }

    // Compute heights (distance from leaves) for all nodes
    // Assumes ultrametric-ish: height = max distance to any leaf
    let heights = compute_heights(tree, root)?;

    // Priority queue of (height, node_id)
    // We want to split the node with the largest height
    use std::cmp::Ordering;
    use std::collections::BinaryHeap;

    #[derive(PartialEq)]
    struct NodeHeight {
        h: f64,
        id: NodeId,
    }

    impl Eq for NodeHeight {}

    impl PartialOrd for NodeHeight {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }

    impl Ord for NodeHeight {
        fn cmp(&self, other: &Self) -> Ordering {
            self.h.partial_cmp(&other.h).unwrap_or(Ordering::Equal)
        }
    }

    let mut leaves_clusters = Vec::new();
    let mut active_queue = BinaryHeap::new();
    active_queue.push(NodeHeight {
        h: heights[&root],
        id: root,
    });

    // Total clusters = leaves_clusters.len() + active_queue.len()
    // Initially 0 + 1 = 1.

    while leaves_clusters.len() + active_queue.len() < k {
        // Try to pick the best candidate to split
        if let Some(top) = active_queue.pop() {
            let node = tree
                .get_node(top.id)
                .ok_or_else(|| anyhow::anyhow!("node {} not found", top.id))?;
            if node.children.is_empty() {
                // Cannot split a leaf, put it back into finished list
                leaves_clusters.push(top);
            } else {
                // Split this node: Add children to queue
                for &child in &node.children {
                    active_queue.push(NodeHeight {
                        h: heights[&child],
                        id: child,
                    });
                }
            }
        } else {
            // No more splittable nodes
            break;
        }
    }

    let mut final_roots = Vec::new();
    for item in leaves_clusters {
        final_roots.push(item.id);
    }
    for item in active_queue {
        final_roots.push(item.id);
    }

    assign_clusters(tree, final_roots)
}

/// Cut tree at specific height (distance from leaves)
pub fn cut_height(tree: &Tree, h: f64) -> Result<Partition> {
    let root = tree
        .get_root()
        .ok_or_else(|| anyhow::anyhow!("Tree has no root"))?;
    let heights = compute_heights(tree, root)?;

    let mut clusters = Vec::new();
    let mut stack = vec![root];

    while let Some(node_id) = stack.pop() {
        let height = heights[&node_id];
        let node = tree
            .get_node(node_id)
            .ok_or_else(|| anyhow::anyhow!("node {} not found", node_id))?;

        if height <= h {
            // This node is below threshold, it forms a cluster
            clusters.push(node_id);
        } else {
            // Node is above threshold
            if node.children.is_empty() {
                // Leaf but height > h? Should not happen if h >= 0, as leaf height is 0.
                clusters.push(node_id);
            } else {
                // Continue down
                for &child in &node.children {
                    stack.push(child);
                }
            }
        }
    }

    assign_clusters(tree, clusters)
}

/// Cut tree at specific distance from root
pub fn cut_root_dist(tree: &Tree, d: f64) -> Result<Partition> {
    let root = tree
        .get_root()
        .ok_or_else(|| anyhow::anyhow!("Tree has no root"))?;

    // We traverse from root.
    // If current_dist + edge_len >= d, then the edge crosses the threshold.
    // The child node represents the cluster (or rather, the subtree at child).
    // Note: If root itself is already beyond d (unlikely if d>0), then root is cluster?

    let mut clusters = Vec::new();
    let mut stack = vec![(root, 0.0)]; // id, current_dist

    while let Some((node_id, dist)) = stack.pop() {
        // If we are already past distance, this shouldn't happen with the logic below,
        // unless root starts past d.
        if dist >= d {
            clusters.push(node_id);
            continue;
        }

        let node = tree
            .get_node(node_id)
            .ok_or_else(|| anyhow::anyhow!("node {} not found", node_id))?;

        if node.children.is_empty() {
            // Leaf reached before cut distance. It's a cluster.
            clusters.push(node_id);
        } else {
            for &child in &node.children {
                let child_node = tree
                    .get_node(child)
                    .ok_or_else(|| anyhow::anyhow!("node {} not found", child))?;
                let len = child_node.length.unwrap_or(0.0);
                let child_dist = dist + len;

                if child_dist >= d {
                    // The edge to child crosses the threshold (or lands exactly on it)
                    // So 'child' becomes the root of a new cluster
                    clusters.push(child);
                } else {
                    // Still within distance, continue traversing
                    stack.push((child, child_dist));
                }
            }
        }
    }

    assign_clusters(tree, clusters)
}

#[cfg(test)]
mod tests {
    use super::*;
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
    fn test_cut_k_two_clusters() {
        // Tree: ((A:1,B:1):1,C:1);
        // K=2 should split at root, giving {A,B} and {C}
        let tree = parse_tree("((A:1,B:1):1,C:1);");
        let part = cut_k(&tree, 2).unwrap();
        let clusters = cluster_names(&part, &tree);
        assert_eq!(clusters, vec![vec!["A", "B"], vec!["C"]]);
    }

    #[test]
    fn test_cut_k_all_singletons() {
        let tree = parse_tree("((A:1,B:1):1,C:1);");
        let part = cut_k(&tree, 3).unwrap();
        let mut clusters = cluster_names(&part, &tree);
        clusters.sort();
        assert_eq!(clusters, vec![vec!["A"], vec!["B"], vec!["C"]]);
    }

    #[test]
    fn test_cut_height_all_singletons() {
        let tree = parse_tree("((A:1,B:1):1,C:1);");
        // h=0.5 cuts all internal nodes
        let part = cut_height(&tree, 0.5).unwrap();
        let mut clusters = cluster_names(&part, &tree);
        clusters.sort();
        assert_eq!(clusters, vec![vec!["A"], vec!["B"], vec!["C"]]);
    }

    #[test]
    fn test_cut_height_one_cluster() {
        let tree = parse_tree("((A:1,B:1):1,C:1);");
        // h=2.0 keeps root (height 2)
        let part = cut_height(&tree, 2.0).unwrap();
        let clusters = cluster_names(&part, &tree);
        assert_eq!(clusters, vec![vec!["A", "B", "C"]]);
    }

    #[test]
    fn test_cut_root_dist_two_clusters() {
        let tree = parse_tree("((A:1,B:1):1,C:1);");
        // d=0.5 cuts both edges from the root, so D=(A,B) and C become cluster roots.
        let part = cut_root_dist(&tree, 0.5).unwrap();
        let mut clusters = cluster_names(&part, &tree);
        clusters.sort();
        assert_eq!(clusters, vec![vec!["A", "B"], vec!["C"]]);
    }

    #[test]
    fn test_cut_root_dist_all_singletons() {
        let tree = parse_tree("((A:1,B:1):1,C:1);");
        // d=1.5 is beyond the A/B internal edges, so every leaf becomes its own cluster.
        let part = cut_root_dist(&tree, 1.5).unwrap();
        let mut clusters = cluster_names(&part, &tree);
        clusters.sort();
        assert_eq!(clusters, vec![vec!["A"], vec!["B"], vec!["C"]]);
    }

    #[test]
    fn test_cut_k_zero_rejected() {
        let tree = parse_tree("(A,B);");
        assert!(cut_k(&tree, 0).is_err());
    }
}
