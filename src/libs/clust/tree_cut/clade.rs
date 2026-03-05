use super::{assign_clusters, Partition};
use crate::libs::phylo::node::NodeId;
use crate::libs::phylo::tree::Tree;
use std::collections::HashMap;

/// TreeCluster: Max pairwise distance in clade <= threshold
pub fn cut_max_clade(tree: &Tree, threshold: f64) -> Result<Partition, String> {
    let root = tree.get_root().ok_or("Tree has no root")?;

    // We need to compute diameters Bottom-Up.
    // Diameter of T(u) is max(Diameter(T(v)), Diameter(T(w)), MaxPathPassingThrough(u))
    // MaxPathPassingThrough(u) = Sum of two largest (MaxDepth(child) + len(u->child)).

    let mut cluster_roots = Vec::new();
    let mut diameters = HashMap::new();
    let mut max_depths = HashMap::new(); // Max distance from u to descendant leaf

    // Post-order for calculation
    let post_order = crate::libs::phylo::tree::traversal::postorder(tree, root);

    for id in post_order {
        let node = tree.get_node(id).unwrap();
        if node.children.is_empty() {
            max_depths.insert(id, 0.0);
            diameters.insert(id, 0.0);
        } else {
            let mut depths = Vec::new();
            let mut child_diams = Vec::new();

            for &child in &node.children {
                let child_node = tree.get_node(child).unwrap();
                let len = child_node.length.unwrap_or(0.0);
                let d = max_depths[&child];
                depths.push(d + len);
                child_diams.push(diameters[&child]);
            }

            // Max depth
            let my_max_depth = depths.iter().cloned().fold(0.0, f64::max);
            max_depths.insert(id, my_max_depth);

            // Diameter
            // 1. Max of children diameters
            let max_child_diam = child_diams.iter().cloned().fold(0.0, f64::max);

            // 2. Path through u: Sum of two largest depths
            let mut sorted_depths = depths;
            sorted_depths.sort_by(|a, b| b.partial_cmp(a).unwrap());
            let path_thru_u = if sorted_depths.len() >= 2 {
                sorted_depths[0] + sorted_depths[1]
            } else {
                0.0
            };

            let my_diam = max_child_diam.max(path_thru_u);
            diameters.insert(id, my_diam);
        }
    }

    // Top-Down to pick maximal clusters
    let mut stack = vec![root];
    while let Some(id) = stack.pop() {
        let d = diameters[&id];
        if d <= threshold {
            cluster_roots.push(id);
        } else {
            let node = tree.get_node(id).unwrap();
            for &child in &node.children {
                stack.push(child);
            }
        }
    }

    assign_clusters(tree, cluster_roots)
}

/// Cut tree where average pairwise distance in cluster <= threshold.
pub fn cut_avg_clade(tree: &Tree, threshold: f64) -> Result<Partition, String> {
    let root = tree.get_root().ok_or("Tree has no root")?;

    // 1. Compute avg distances
    let avg_dists = crate::libs::phylo::tree::stat::compute_avg_clade_distances(tree);

    // 2. Top-down greedy cut
    let mut clusters = Vec::new();
    let mut stack = vec![root];

    while let Some(node_id) = stack.pop() {
        let avg_dist = *avg_dists.get(&node_id).unwrap_or(&0.0);

        if avg_dist <= threshold {
            clusters.push(node_id);
        } else {
            // Split
            if let Some(node) = tree.get_node(node_id) {
                if node.children.is_empty() {
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

/// Cut tree based on median pairwise distance in clade.
pub fn cut_med_clade(tree: &Tree, threshold: f64) -> Result<Partition, String> {
    let root = tree.get_root().ok_or("Tree has no root")?;
    let mut clusters = Vec::new();

    // We will simulate the bottom-up pass using post-order traversal.
    let mut post_order = Vec::new();
    let mut visit_stack = vec![root];
    while let Some(u) = visit_stack.pop() {
        post_order.push(u);
        if let Some(node) = tree.get_node(u) {
            for &child in &node.children {
                visit_stack.push(child);
            }
        }
    }
    // visit_stack pop order is Pre-Order (Parent, Child).
    // Reversed, it becomes Post-Order (Child, Parent).
    let post_order: Vec<NodeId> = post_order.into_iter().rev().collect();

    // Map NodeId -> (LeafDists, PairDists)
    // LeafDists: sorted list of distances from leaves to current node
    // PairDists: sorted list of all pairwise distances in subtree
    // Using simple Vec and sorting for now. Merging sorted lists is better but more code.
    struct NodeStat {
        leaf_dists: Vec<f64>,
        pair_dists: Vec<f64>,
        median: f64,
    }

    let mut stats: HashMap<NodeId, NodeStat> = HashMap::new();

    for &u in &post_order {
        let node = tree.get_node(u).unwrap();
        if node.children.is_empty() {
            // Leaf
            stats.insert(
                u,
                NodeStat {
                    leaf_dists: vec![0.0],
                    pair_dists: Vec::new(),
                    median: 0.0,
                },
            );
        } else {
            // Merge children
            let mut my_leaf_dists = Vec::new();
            let mut my_pair_dists = Vec::new();

            // 1. Update leaf distances and accumulate pair distances from subtrees
            let mut child_leaf_dists_list = Vec::new();

            for &v in &node.children {
                let child_stat = stats.get(&v).unwrap();
                let len = tree.get_node(v).unwrap().length.unwrap_or(0.0);

                // Add child's pair dists
                my_pair_dists.extend_from_slice(&child_stat.pair_dists);

                // Shift child's leaf dists by edge length
                let shifted_leaf_dists: Vec<f64> =
                    child_stat.leaf_dists.iter().map(|d| d + len).collect();
                child_leaf_dists_list.push(shifted_leaf_dists.clone());
                my_leaf_dists.extend(shifted_leaf_dists);
            }

            // 2. Compute cross distances between subtrees
            for i in 0..child_leaf_dists_list.len() {
                for j in (i + 1)..child_leaf_dists_list.len() {
                    let list_a = &child_leaf_dists_list[i];
                    let list_b = &child_leaf_dists_list[j];

                    // Cross product: O(La * Lb)
                    for &da in list_a {
                        for &db in list_b {
                            my_pair_dists.push(da + db);
                        }
                    }
                }
            }

            // 3. Compute median
            if my_pair_dists.is_empty() {
                stats.insert(
                    u,
                    NodeStat {
                        leaf_dists: my_leaf_dists,
                        pair_dists: my_pair_dists,
                        median: 0.0,
                    },
                );
            } else {
                // Sort to find median
                my_pair_dists.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

                let mid = my_pair_dists.len() / 2;
                let median = if my_pair_dists.len() % 2 == 1 {
                    my_pair_dists[mid]
                } else {
                    (my_pair_dists[mid - 1] + my_pair_dists[mid]) / 2.0
                };

                #[cfg(test)]
                {
                    // Get node name for debug
                    let name = tree.get_node(u).unwrap().name.as_deref().unwrap_or("?");
                    println!(
                        "DEBUG: Node {} ({}): Median {}, PairDists {:?}",
                        u, name, median, my_pair_dists
                    );
                }

                stats.insert(
                    u,
                    NodeStat {
                        leaf_dists: my_leaf_dists,
                        pair_dists: my_pair_dists,
                        median,
                    },
                );
            }
        }
    }

    // Top-down selection
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(root);

    while let Some(u) = queue.pop_front() {
        let stat = stats.get(&u).unwrap();
        // Use epsilon for float comparison? TreeCluster uses <=
        if stat.median <= threshold + 1e-9 {
            clusters.push(u);
        } else {
            if let Some(node) = tree.get_node(u) {
                for &child in &node.children {
                    queue.push_back(child);
                }
            }
        }
    }

    assign_clusters(tree, clusters)
}

/// Cut tree based on sum of branch lengths in clade.
pub fn cut_sum_branch(tree: &Tree, threshold: f64) -> Result<Partition, String> {
    let root = tree.get_root().ok_or("Tree has no root")?;
    let mut clusters = Vec::new();

    // Map NodeId -> Subtree Sum Branch Length
    let mut sums: HashMap<NodeId, f64> = HashMap::new();

    // Post-order traversal for bottom-up calculation
    let mut post_order = Vec::new();
    let mut visit_stack = vec![root];
    while let Some(u) = visit_stack.pop() {
        post_order.push(u);
        if let Some(node) = tree.get_node(u) {
            for &child in &node.children {
                visit_stack.push(child);
            }
        }
    }
    let post_order: Vec<NodeId> = post_order.into_iter().rev().collect();

    for &u in &post_order {
        let node = tree.get_node(u).unwrap();
        if node.children.is_empty() {
            // Leaf has 0 internal branch length
            sums.insert(u, 0.0);
        } else {
            let mut sum = 0.0;
            for &v in &node.children {
                let child_sum = sums.get(&v).unwrap();
                let len = tree.get_node(v).unwrap().length.unwrap_or(0.0);
                sum += child_sum + len;
            }
            sums.insert(u, sum);
        }
    }

    // Top-down selection
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(root);

    while let Some(u) = queue.pop_front() {
        let sum = sums.get(&u).unwrap();
        // Use epsilon?
        if *sum <= threshold + 1e-9 {
            clusters.push(u);
        } else {
            if let Some(node) = tree.get_node(u) {
                for &child in &node.children {
                    queue.push_back(child);
                }
            }
        }
    }

    assign_clusters(tree, clusters)
}
