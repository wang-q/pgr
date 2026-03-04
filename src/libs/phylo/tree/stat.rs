use super::Tree;
use crate::libs::phylo::node::NodeId;
use std::collections::{BTreeMap, HashMap, VecDeque};

/// Get IDs of all leaves in subtree rooted at `id`.
pub fn get_leaves(tree: &Tree, id: NodeId) -> Vec<NodeId> {
    let mut leaves = Vec::new();
    let mut stack = vec![id];

    while let Some(curr) = stack.pop() {
        if let Some(node) = tree.get_node(curr) {
            if node.children.is_empty() {
                leaves.push(curr);
            } else {
                for &child in &node.children {
                    stack.push(child);
                }
            }
        }
    }
    leaves
}

/// Get names of all leaves in subtree.
pub fn get_leaf_names(tree: &Tree, id: NodeId) -> Vec<Option<String>> {
    get_leaves(tree, id)
        .into_iter()
        .map(|leaf_id| tree.get_node(leaf_id).and_then(|n| n.name.clone()))
        .collect()
}

/// Check if tree is binary (all internal nodes have degree 2).
/// Note: Root can have degree 2 (bifurcating) or 3 (unrooted representation) or more.
/// This checks if *children count* is 2 for all internal nodes.
pub fn is_binary(tree: &Tree) -> bool {
    tree.nodes
        .iter()
        .all(|n| n.deleted || n.children.is_empty() || n.children.len() == 2)
}

/// Check if the tree is rooted (root node has degree 2).
/// Unrooted trees typically have a trifurcating root (degree >= 3).
pub fn is_rooted(tree: &Tree) -> bool {
    if let Some(root_id) = tree.get_root() {
        if let Some(node) = tree.get_node(root_id) {
            return node.children.len() == 2;
        }
    }
    false
}

/// Calculate diameter (longest path between any two nodes).
pub fn diameter(tree: &Tree, weighted: bool) -> f64 {
    // 2-pass BFS/DFS to find diameter.
    // 1. Find furthest node from Root (A)
    // 2. Find furthest node from A (B)
    // Dist(A, B) is diameter.

    let root = match tree.get_root() {
        Some(r) => r,
        None => return 0.0,
    };

    // Helper: BFS to find (node, distance)
    let get_furthest = |start: NodeId| -> (NodeId, f64) {
        let mut max_dist = 0.0;
        let mut furthest_node = start;
        let mut visited = HashMap::new();
        let mut queue = std::collections::VecDeque::new();

        visited.insert(start, 0.0);
        queue.push_back(start);

        while let Some(u) = queue.pop_front() {
            let d = *visited.get(&u).unwrap();
            if d > max_dist {
                max_dist = d;
                furthest_node = u;
            }

            // Neighbors: children + parent
            let node = tree.get_node(u).unwrap();
            let mut neighbors = node.children.clone();
            if let Some(p) = node.parent {
                neighbors.push(p);
            }

            for v in neighbors {
                if !visited.contains_key(&v) {
                    // Edge weight
                    let weight = if weighted {
                        // Edge is between u and v. One is parent of other.
                        let v_node = tree.get_node(v).unwrap();
                        let u_node = tree.get_node(u).unwrap();
                        if v_node.parent == Some(u) {
                            v_node.length.unwrap_or(0.0)
                        } else {
                            // u is child of v
                            u_node.length.unwrap_or(0.0)
                        }
                    } else {
                        1.0
                    };

                    visited.insert(v, d + weight);
                    queue.push_back(v);
                }
            }
        }
        (furthest_node, max_dist)
    };

    let (node_a, _) = get_furthest(root);
    let (_, diam) = get_furthest(node_a);
    diam
}

/// Computes the number of cherries in a tree.
/// A cherry is a pair of leaves that share a common parent.
pub fn cherries(tree: &Tree) -> usize {
    let mut count = 0;
    for node in &tree.nodes {
        if !node.deleted && node.children.len() == 2 {
            // Check if both children are leaves
            // We need to handle potential deleted children if that's possible,
            // but tree.nodes iteration includes deleted ones unless checked.
            // node.children contains IDs.

            let is_cherry = node
                .children
                .iter()
                .all(|&child_id| match tree.get_node(child_id) {
                    Some(child) => child.children.is_empty(),
                    None => false,
                });

            if is_cherry {
                count += 1;
            }
        }
    }
    count
}

/// Computes the Sackin index.
/// The Sackin index is the sum of depths of all leaves.
/// Smaller Sackin index means a more balanced tree.
pub fn sackin(tree: &Tree) -> usize {
    let root = match tree.get_root() {
        Some(r) => r,
        None => return 0,
    };

    let mut sum_depth = 0;
    let mut stack = vec![(root, 0)];

    while let Some((node_id, depth)) = stack.pop() {
        if let Some(node) = tree.get_node(node_id) {
            if node.children.is_empty() {
                sum_depth += depth;
            } else {
                for &child in &node.children {
                    stack.push((child, depth + 1));
                }
            }
        }
    }
    sum_depth
}

/// Computes the Colless index.
/// The Colless index is the sum of absolute differences between the number of leaves
/// in the left and right subtrees of each internal node.
/// Only defined for binary trees. Returns None if tree is not binary.
pub fn colless(tree: &Tree) -> Option<usize> {
    if !is_binary(tree) {
        return None;
    }
    let root = tree.get_root()?;

    // Post-order to compute leaf counts
    let post_order = crate::libs::phylo::tree::traversal::postorder(tree, root);
    let mut leaf_counts = HashMap::new();
    let mut index = 0;

    for id in post_order {
        if let Some(node) = tree.get_node(id) {
            if node.children.is_empty() {
                leaf_counts.insert(id, 1);
            } else {
                let mut count = 0;
                let mut children_counts = Vec::new();
                for &child in &node.children {
                    let c = *leaf_counts.get(&child).unwrap_or(&0);
                    count += c;
                    children_counts.push(c);
                }
                leaf_counts.insert(id, count);

                if children_counts.len() == 2 {
                    let diff = (children_counts[0] as isize - children_counts[1] as isize).abs();
                    index += diff as usize;
                }
            }
        }
    }
    Some(index)
}

/// Get all node names.
pub fn get_names(tree: &Tree) -> Vec<String> {
    tree.nodes
        .iter()
        .filter(|n| !n.deleted)
        .filter_map(|n| n.name.clone())
        .collect()
}

/// Get mapping from name to NodeId.
pub fn get_name_id(tree: &Tree) -> BTreeMap<String, NodeId> {
    let mut map = BTreeMap::new();
    for node in &tree.nodes {
        if !node.deleted {
            if let Some(name) = &node.name {
                map.insert(name.clone(), node.id);
            }
        }
    }
    map
}

/// Get values for a specific property key for all nodes.
pub fn get_property_values(tree: &Tree, key: &str) -> BTreeMap<NodeId, String> {
    tree.nodes
        .iter()
        .filter(|n| !n.deleted)
        .filter_map(|n| {
            n.properties
                .as_ref()
                .and_then(|p| p.get(key))
                .map(|val| (n.id, val.clone()))
        })
        .collect()
}

/// Find the node with the longest edge length.
pub fn get_node_with_longest_edge(tree: &Tree) -> Option<NodeId> {
    let mut max_len = f64::NEG_INFINITY;
    let mut max_node = None;

    for node in &tree.nodes {
        if !node.deleted {
            if let Some(len) = node.length {
                if len > max_len {
                    max_len = len;
                    max_node = Some(node.id);
                }
            }
        }
    }
    max_node
}

/// Compute node heights (distance from leaves).
/// Assumes ultrametric-like tree (takes max distance to leaves).
pub fn compute_node_heights(tree: &Tree) -> HashMap<NodeId, f64> {
    let mut heights = HashMap::new();
    if let Some(root) = tree.get_root() {
        let ids = super::traversal::postorder(tree, root);

        for id in ids {
            if let Some(node) = tree.get_node(id) {
                if node.children.is_empty() {
                    heights.insert(id, 0.0);
                } else {
                    let mut max_h = 0.0;
                    for &child in &node.children {
                        let child_h = *heights.get(&child).unwrap_or(&0.0);
                        let len = tree.get_node(child).and_then(|n| n.length).unwrap_or(0.0);
                        let h = child_h + len;
                        if h > max_h {
                            max_h = h;
                        }
                    }
                    heights.insert(id, max_h);
                }
            }
        }
    }
    heights
}

/// Calculate inconsistent coefficients for all internal nodes.
/// Returns a map from NodeId to coefficient.
///
/// Coefficient = (height(u) - mean(heights in subtree)) / std(heights in subtree)
/// Subtree includes u and descendants up to depth `d`.
pub fn calculate_inconsistency(
    tree: &Tree,
    node_heights: &HashMap<NodeId, f64>,
    depth: usize,
) -> HashMap<NodeId, f64> {
    let mut coeffs = HashMap::new();

    // We only care about internal nodes
    for (&id, &h) in node_heights {
        if let Some(node) = tree.get_node(id) {
            if node.children.is_empty() {
                continue;
            }

            // Collect heights in subtree up to depth
            let mut collected = Vec::new();
            let mut queue = VecDeque::new();
            queue.push_back((id, 0)); // (node_id, current_depth)

            while let Some((curr, d)) = queue.pop_front() {
                // Only include internal nodes (links)
                let is_internal = tree
                    .get_node(curr)
                    .map(|n| !n.children.is_empty())
                    .unwrap_or(false);

                if is_internal {
                    if let Some(&curr_h) = node_heights.get(&curr) {
                        collected.push(curr_h);
                    }
                }

                if d < depth {
                    if let Some(curr_node) = tree.get_node(curr) {
                        for &child in &curr_node.children {
                            let child_internal = tree
                                .get_node(child)
                                .map(|n| !n.children.is_empty())
                                .unwrap_or(false);
                            if child_internal {
                                queue.push_back((child, d + 1));
                            }
                        }
                    }
                }
            }

            if collected.is_empty() {
                coeffs.insert(id, 0.0);
                continue;
            }

            let n = collected.len() as f64;
            let mean = collected.iter().sum::<f64>() / n;

            if n <= 1.0 {
                coeffs.insert(id, 0.0);
                continue;
            }

            // Use sample variance (ddof=1) to match likely SciPy behavior (which passes 0.8 threshold for n=2)
            let variance = collected.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0);

            let std = if variance > 1e-9 {
                variance.sqrt()
            } else {
                0.0
            };

            if std == 0.0 {
                coeffs.insert(id, 0.0);
            } else {
                coeffs.insert(id, (h - mean) / std);
            }
        }
    }
    coeffs
}

/// Compute average pairwise distance between leaves for each cluster (subtree).
/// Returns a map from NodeId to average distance.
pub fn compute_avg_clade_distances(tree: &Tree) -> HashMap<NodeId, f64> {
    let mut avg_dists = HashMap::new();
    // Helper stats: (num_leaves, sum_leaf_dist_to_node, sum_pair_dist)
    let mut stats: HashMap<NodeId, (usize, f64, f64)> = HashMap::new();

    if let Some(root) = tree.get_root() {
        let post_order = super::traversal::postorder(tree, root);
        for id in post_order {
            if let Some(node) = tree.get_node(id) {
                if node.children.is_empty() {
                    stats.insert(id, (1, 0.0, 0.0));
                    avg_dists.insert(id, 0.0);
                } else {
                    let mut num_leaves = 0;
                    let mut sum_leaf_dist = 0.0;
                    let mut sum_pair_dist = 0.0;

                    let mut child_stats = Vec::new();

                    for &child in &node.children {
                        let (c_n, c_leaf_dist, c_pair_dist) =
                            *stats.get(&child).unwrap_or(&(0, 0.0, 0.0));
                        let len = tree.get_node(child).and_then(|n| n.length).unwrap_or(0.0);
                        let c_leaf_dist_to_u = c_leaf_dist + (c_n as f64) * len;

                        child_stats.push((c_n, c_leaf_dist_to_u));

                        num_leaves += c_n;
                        sum_leaf_dist += c_leaf_dist_to_u;
                        sum_pair_dist += c_pair_dist;
                    }

                    // Add cross-cluster pairs
                    for (c_n, c_leaf_dist_to_u) in child_stats {
                        sum_pair_dist += c_leaf_dist_to_u * (num_leaves - c_n) as f64;
                    }

                    stats.insert(id, (num_leaves, sum_leaf_dist, sum_pair_dist));

                    if num_leaves > 1 {
                        let pairs = (num_leaves * (num_leaves - 1)) as f64 / 2.0;
                        avg_dists.insert(id, sum_pair_dist / pairs);
                    } else {
                        avg_dists.insert(id, 0.0);
                    }
                }
            }
        }
    }
    avg_dists
}
