//! Tree partitioning algorithms.
//!
//! This module provides methods to cut a phylogenetic tree into clusters based on various criteria.
//! The goal is to partition the leaves into disjoint sets (clusters).
//!
//! # Supported Methods
//!
//! *   **K**: Cut the tree into exactly K clusters.
//!     This is a top-down approach that iteratively splits the cluster with the largest height
//!     until K clusters are formed.
//!
//! *   **Height**: Cut the tree at a specific distance from the leaves (ultrametric assumption).
//!     Any node with height > threshold is split; nodes with height <= threshold form clusters.
//!
//! *   **RootDist**: Cut the tree at a specific distance from the root.
//!     This simulates a "timeline" cut. Branches crossing the threshold are cut, and the
//!     subtrees below become clusters.
//!
//! *   **MaxClade**: Ensure the maximum pairwise distance (diameter) within any cluster
//!     does not exceed a threshold. This corresponds to the `MaxClade` criterion in TreeCluster.
//!

use super::Tree;
use crate::libs::phylo::node::NodeId;
use std::collections::HashMap;

/// Cut strategies for tree partitioning.
#[derive(Debug, Clone, Copy)]
pub enum Method {
    /// Cut into exactly K clusters.
    ///
    /// Iteratively splits the cluster with the largest height (distance to farthest leaf).
    K(usize),

    /// Cut at specific height (distance from leaves).
    ///
    /// Useful for ultrametric trees where height represents time/divergence.
    Height(f64),

    /// Cut at specific distance from root.
    ///
    /// Useful for defining clusters based on divergence from a common ancestor (root).
    RootDist(f64),

    /// TreeCluster: Max pairwise distance in clade <= threshold.
    ///
    /// Ensures that for every cluster, the maximum distance between any two leaves
    /// in that cluster is at most `threshold`.
    MaxClade(f64),

    /// SciPy: Inconsistent coefficient <= threshold.
    ///
    /// Splits nodes if their inconsistency coefficient > threshold.
    /// Requires checking inconsistency of all descendants.
    /// Parameters: (threshold, depth).
    Inconsistent(f64, usize),
}

/// Result of a cut operation.
pub struct Partition {
    /// Map from Leaf NodeId to Cluster ID (1-based).
    pub assignment: HashMap<NodeId, usize>,
    /// Total number of clusters formed.
    pub num_clusters: usize,
}

impl Partition {
    /// Create a new empty partition.
    pub fn new() -> Self {
        Self {
            assignment: HashMap::new(),
            num_clusters: 0,
        }
    }

    /// Get members of each cluster.
    ///
    /// Returns a map where keys are Cluster IDs (1-based) and values are lists of Leaf NodeIds.
    pub fn get_clusters(&self) -> HashMap<usize, Vec<NodeId>> {
        let mut clusters = HashMap::new();
        for (&node_id, &cluster_id) in &self.assignment {
            clusters
                .entry(cluster_id)
                .or_insert_with(Vec::new)
                .push(node_id);
        }
        clusters
    }
}

impl Default for Partition {
    fn default() -> Self {
        Self::new()
    }
}

/// Cut the tree according to the specified method.
///
/// # Arguments
///
/// *   `tree` - The phylogenetic tree to cut.
/// *   `method` - The cutting method/strategy.
///
/// # Returns
///
/// A `Result` containing the `Partition` or an error message.
///
/// # Examples
///
/// ```ignore
/// use pgr::libs::phylo::tree::cut::{cut, Method};
///
/// let tree = ...; // Load or create tree
/// let partition = cut(&tree, Method::K(3)).unwrap();
/// println!("Formed {} clusters", partition.num_clusters);
/// ```
pub fn cut(tree: &Tree, method: Method) -> Result<Partition, String> {
    if tree.is_empty() {
        return Ok(Partition::new());
    }

    match method {
        Method::K(k) => cut_k(tree, k),
        Method::Height(h) => cut_height(tree, h),
        Method::RootDist(d) => cut_root_dist(tree, d),
        Method::MaxClade(t) => cut_max_clade(tree, t),
        Method::Inconsistent(t, d) => cut_inconsistent(tree, t, d),
    }
}

/// Cut tree into K clusters
fn cut_k(tree: &Tree, k: usize) -> Result<Partition, String> {
    if k == 0 {
        return Err("K must be >= 1".to_string());
    }

    let root = tree.get_root().ok_or("Tree has no root")?;
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
            self.h.partial_cmp(&other.h)
        }
    }

    impl Ord for NodeHeight {
        fn cmp(&self, other: &Self) -> Ordering {
            self.partial_cmp(other).unwrap_or(Ordering::Equal)
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
            let node = tree.get_node(top.id).unwrap();
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
fn cut_height(tree: &Tree, h: f64) -> Result<Partition, String> {
    let root = tree.get_root().ok_or("Tree has no root")?;
    let heights = compute_heights(tree, root)?;

    let mut clusters = Vec::new();
    let mut stack = vec![root];

    while let Some(node_id) = stack.pop() {
        let height = heights[&node_id];
        let node = tree.get_node(node_id).unwrap();

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
fn cut_root_dist(tree: &Tree, d: f64) -> Result<Partition, String> {
    let root = tree.get_root().ok_or("Tree has no root")?;

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

        let node = tree.get_node(node_id).unwrap();

        if node.children.is_empty() {
            // Leaf reached before cut distance. It's a cluster.
            clusters.push(node_id);
        } else {
            for &child in &node.children {
                let child_node = tree.get_node(child).unwrap();
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

/// TreeCluster: Max pairwise distance in clade <= threshold
fn cut_max_clade(tree: &Tree, threshold: f64) -> Result<Partition, String> {
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

/// Cut tree based on inconsistent coefficient threshold.
///
/// Corresponds to SciPy `fcluster(criterion='inconsistent')`.
/// A node `u` forms a cluster if `inconsistent(v) <= threshold` for all `v` in subtree of `u` (including `u`).
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

// --- Helpers ---

/// Compute max distance from each node to its leaves
fn compute_heights(tree: &Tree, root: NodeId) -> Result<HashMap<NodeId, f64>, String> {
    let mut heights = HashMap::new();
    let post_order = crate::libs::phylo::tree::traversal::postorder(tree, root);

    for id in post_order {
        let node = tree.get_node(id).unwrap();
        if node.children.is_empty() {
            heights.insert(id, 0.0);
        } else {
            let mut max_h = 0.0;
            for &child in &node.children {
                let child_node = tree.get_node(child).unwrap();
                let len = child_node.length.unwrap_or(0.0); // If None, assume 0
                let h = heights[&child] + len;
                if h > max_h {
                    max_h = h;
                }
            }
            heights.insert(id, max_h);
        }
    }
    Ok(heights)
}

/// Assign leaves to clusters based on cluster roots
fn assign_clusters(tree: &Tree, cluster_roots: Vec<NodeId>) -> Result<Partition, String> {
    let mut part = Partition::new();
    let mut cluster_id = 0;

    for root in cluster_roots {
        cluster_id += 1;
        let leaves = crate::libs::phylo::tree::stat::get_leaves(tree, root);
        for leaf in leaves {
            part.assignment.insert(leaf, cluster_id);
        }
    }

    part.num_clusters = cluster_id;
    Ok(part)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::phylo::tree::Tree;

    fn create_test_tree() -> Tree {
        // ((A:1,B:1)D:1,C:1)E;
        // Leaf heights: A=0, B=0, C=0
        // Node D: height=1
        // Root E: height=2
        // Distances from root E: D=1, C=1, A=2, B=2
        // MaxClade (Diameter): D=2, E=3
        let mut tree = Tree::new();
        let a = tree.add_node();
        tree.get_node_mut(a).unwrap().name = Some("A".to_string());
        let b = tree.add_node();
        tree.get_node_mut(b).unwrap().name = Some("B".to_string());
        let c = tree.add_node();
        tree.get_node_mut(c).unwrap().name = Some("C".to_string());
        let d = tree.add_node();
        tree.get_node_mut(d).unwrap().name = Some("D".to_string());
        let e = tree.add_node();
        tree.get_node_mut(e).unwrap().name = Some("E".to_string());

        tree.set_root(e);

        // E -> D, C
        tree.add_child(e, d).unwrap();
        tree.add_child(e, c).unwrap();
        tree.get_node_mut(d).unwrap().length = Some(1.0);
        tree.get_node_mut(c).unwrap().length = Some(1.0);

        // D -> A, B
        tree.add_child(d, a).unwrap();
        tree.add_child(d, b).unwrap();
        tree.get_node_mut(a).unwrap().length = Some(1.0);
        tree.get_node_mut(b).unwrap().length = Some(1.0);

        tree
    }

    #[test]
    fn test_cut_k() {
        let tree = create_test_tree();

        // K=1 -> {A,B,C} (1 cluster)
        let p1 = cut(&tree, Method::K(1)).unwrap();
        assert_eq!(p1.num_clusters, 1);

        // K=2 -> {A,B} (from D), {C} (from C)
        let p2 = cut(&tree, Method::K(2)).unwrap();
        assert_eq!(p2.num_clusters, 2);
        let clusters = p2.get_clusters();
        assert_eq!(clusters.len(), 2);
        // We don't know cluster IDs, but we know grouping
        let a = 0;
        let b = 1;
        let c = 2;
        assert_eq!(p2.assignment[&a], p2.assignment[&b]);
        assert_ne!(p2.assignment[&a], p2.assignment[&c]);

        // K=3 -> {A}, {B}, {C}
        let p3 = cut(&tree, Method::K(3)).unwrap();
        assert_eq!(p3.num_clusters, 3);
        assert_ne!(p3.assignment[&a], p3.assignment[&b]);
        assert_ne!(p3.assignment[&a], p3.assignment[&c]);
    }

    #[test]
    fn test_cut_height() {
        let tree = create_test_tree();

        // Height=1.5
        // E(2.0) > 1.5 -> split
        // D(1.0) <= 1.5 -> cluster
        // C(0.0) <= 1.5 -> cluster
        // Result: {A,B}, {C}
        let p = cut(&tree, Method::Height(1.5)).unwrap();
        assert_eq!(p.num_clusters, 2);
        let a = 0;
        let b = 1;
        let c = 2;
        assert_eq!(p.assignment[&a], p.assignment[&b]);
        assert_ne!(p.assignment[&a], p.assignment[&c]);

        // Height=0.5
        // D(1.0) > 0.5 -> split
        // A, B, C <= 0.5 -> clusters
        let p2 = cut(&tree, Method::Height(0.5)).unwrap();
        assert_eq!(p2.num_clusters, 3);
    }

    #[test]
    fn test_cut_root_dist() {
        let tree = create_test_tree();

        // Cut at 0.5 from root
        // E -> D (len 1), E -> C (len 1)
        // Both edges cross 0.5
        // Clusters: {A,B}, {C}
        let p = cut(&tree, Method::RootDist(0.5)).unwrap();
        assert_eq!(p.num_clusters, 2);
        let a = 0;
        let b = 1;
        let c = 2;
        assert_eq!(p.assignment[&a], p.assignment[&b]);
        assert_ne!(p.assignment[&a], p.assignment[&c]);

        // Cut at 1.5
        // E->D (1.0) < 1.5 -> traverse D
        // D->A (1.0+1.0=2.0) > 1.5 -> split A
        // D->B (2.0) > 1.5 -> split B
        // E->C (1.0) < 1.5 -> traverse C (leaf) -> C is cluster
        // A, B, C are clusters.
        let p2 = cut(&tree, Method::RootDist(1.5)).unwrap();
        assert_eq!(p2.num_clusters, 3);
    }

    #[test]
    fn test_cut_max_clade() {
        let tree = create_test_tree();

        // Threshold 2.5
        // D diameter = 2.0 <= 2.5 -> D is cluster
        // E diameter = 3.0 > 2.5 -> E split
        // C diameter = 0 <= 2.5 -> C is cluster
        // Result: {A,B}, {C}
        let p = cut(&tree, Method::MaxClade(2.5)).unwrap();
        assert_eq!(p.num_clusters, 2);
        let a = 0;
        let b = 1;
        let c = 2;
        assert_eq!(p.assignment[&a], p.assignment[&b]);
        assert_ne!(p.assignment[&a], p.assignment[&c]);

        // Threshold 1.5
        // D diameter = 2.0 > 1.5 -> split
        // Result: {A}, {B}, {C}
        let p2 = cut(&tree, Method::MaxClade(1.5)).unwrap();
        assert_eq!(p2.num_clusters, 3);
    }
}
