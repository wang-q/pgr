use std::collections::{HashMap, VecDeque};

use super::Tree;
use crate::libs::phylo::node::NodeId;

/// Count the number of cherries in the tree.
///
/// A cherry is an internal node with exactly two leaf children.
pub fn cherries(tree: &Tree) -> usize {
    let mut count = 0;
    for node in &tree.nodes {
        if node.deleted || node.children.is_empty() {
            continue;
        }
        // Cherry: internal node with 2 leaf children
        if node.children.len() == 2 {
            // Children may point to deleted nodes in malformed trees; skip those.
            let Some(c1) = tree.get_node(node.children[0]) else {
                continue;
            };
            let Some(c2) = tree.get_node(node.children[1]) else {
                continue;
            };
            if c1.children.is_empty() && c2.children.is_empty() {
                count += 1;
            }
        }
    }
    count
}

/// Sackin index: sum of depths of all leaves (number of edges from root).
///
/// Returns `None` if the tree has no root.
pub fn sackin(tree: &Tree) -> Option<f64> {
    // Sum of depths of leaves (number of edges from root)
    let root = tree.get_root()?;
    let mut sum_depth = 0.0;
    let mut stack = vec![(root, 0)];

    while let Some((id, depth)) = stack.pop() {
        if let Some(node) = tree.get_node(id) {
            if node.children.is_empty() {
                sum_depth += depth as f64;
            } else {
                for &child in &node.children {
                    stack.push((child, depth + 1));
                }
            }
        }
    }
    Some(sum_depth)
}

/// Colless index: sum of |nL - nR| for all internal nodes.
///
/// Only defined for bifurcating trees. Returns `None` if the tree has no root.
pub fn colless(tree: &Tree) -> Option<f64> {
    // Sum |nL - nR| for all internal nodes. Only defined for bifurcating trees.
    let root = tree.get_root()?;
    let mut leaf_counts: HashMap<NodeId, usize> = HashMap::new();
    let post_order = super::traversal::postorder(tree, root);
    let mut sum_diff = 0.0;

    for &id in &post_order {
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
                    let diff =
                        (children_counts[0] as isize - children_counts[1] as isize).abs() as f64;
                    sum_diff += diff;
                }
            }
        }
    }
    Some(sum_diff)
}

/// Compute average pairwise distance between leaves within each clade (subtree).
///
/// Returns a map from node ID to the average distance for the clade rooted at
/// that node.
pub fn compute_avg_clade_distances(tree: &Tree) -> HashMap<NodeId, f64> {
    struct Stat {
        s: f64,
        d: f64,
        n: usize,
    }

    let mut avg_dists = HashMap::new();
    let mut stats: HashMap<NodeId, Stat> = HashMap::new();

    if let Some(root) = tree.get_root() {
        let post_order = super::traversal::postorder(tree, root);
        for &id in &post_order {
            if let Some(node) = tree.get_node(id) {
                if node.children.is_empty() {
                    stats.insert(
                        id,
                        Stat {
                            s: 0.0,
                            d: 0.0,
                            n: 1,
                        },
                    );
                    avg_dists.insert(id, 0.0);
                } else {
                    let mut s_total = 0.0;
                    let mut d_total = 0.0;
                    let mut n_total = 0;

                    for &child in &node.children {
                        if let Some(st) = stats.get(&child) {
                            let len =
                                super::finite_length(tree.get_node(child).and_then(|n| n.length));
                            let d_child_ext = st.d + st.n as f64 * len;

                            let cross = n_total as f64 * d_child_ext + st.n as f64 * d_total;
                            s_total += st.s + cross;
                            d_total += d_child_ext;
                            n_total += st.n;
                        }
                    }

                    stats.insert(
                        id,
                        Stat {
                            s: s_total,
                            d: d_total,
                            n: n_total,
                        },
                    );

                    if n_total > 1 {
                        let pairs = n_total as f64 * (n_total as f64 - 1.0);
                        avg_dists.insert(id, 2.0 * s_total / pairs);
                    } else {
                        avg_dists.insert(id, 0.0);
                    }
                }
            }
        }
    }
    avg_dists
}

/// Calculate inconsistency scores for internal nodes based on subtree heights.
///
/// For each internal node, examines the heights of its descendants within
/// `depth` levels and computes a z-score of the node's own height relative to
/// its local neighborhood.
pub fn calculate_inconsistency(
    tree: &Tree,
    heights: &HashMap<NodeId, f64>,
    depth: usize,
) -> HashMap<NodeId, f64> {
    let mut inconsistency = HashMap::new();

    // Iterate over all nodes, treat each as root of a small subtree
    // Since we need to access all nodes, we can use the `nodes` vector index if it's dense,
    // but better to traverse.
    // However, the function needs random access to `heights`.

    // We can iterate 0..tree.nodes.len()
    for id in 0..tree.nodes.len() {
        if let Some(node) = tree.get_node(id) {
            if node.children.is_empty() {
                continue;
            }

            let h = *heights.get(&id).unwrap_or(&0.0);
            let mut sub_heights = vec![h];

            let mut queue = VecDeque::new();
            queue.push_back((id, 0));

            while let Some((curr, d)) = queue.pop_front() {
                if d >= depth {
                    continue;
                }
                if let Some(curr_node) = tree.get_node(curr) {
                    for &child in &curr_node.children {
                        if let Some(child_node) = tree.get_node(child) {
                            if !child_node.children.is_empty() {
                                let ch = *heights.get(&child).unwrap_or(&0.0);
                                sub_heights.push(ch);
                                queue.push_back((child, d + 1));
                            }
                        }
                    }
                }
            }

            let n = sub_heights.len();
            let mean = sub_heights.iter().sum::<f64>() / n as f64;
            // Use sample variance (divisor n-1) to match SciPy/standard stats
            let divisor = if n > 1 { (n - 1) as f64 } else { 1.0 };
            let variance = sub_heights.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / divisor;
            let std = variance.sqrt();
            if std == 0.0 {
                inconsistency.insert(id, 0.0);
            } else {
                inconsistency.insert(id, (h - mean) / std);
            }
        }
    }

    inconsistency
}
