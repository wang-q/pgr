use crate::libs::phylo::node::NodeId;
use crate::libs::phylo::tree::Tree;
use std::collections::HashMap;

pub mod clade;
pub mod dynamic;
pub mod inconsistent;
pub mod link;
pub mod method;
pub mod simple;

pub use method::Method;

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

    /// Compute summary statistics for the partition.
    /// Returns (num_clusters, num_singletons, num_non_singletons, max_cluster_size).
    pub fn get_stats(&self) -> (usize, usize, usize, usize) {
        let mut sizes = HashMap::new();
        for &cluster_id in self.assignment.values() {
            *sizes.entry(cluster_id).or_insert(0) += 1;
        }
        let mut singletons = 0;
        let mut non_singletons = 0;
        let mut max_size = 0;
        for &size in sizes.values() {
            if size == 1 {
                singletons += 1;
            } else {
                non_singletons += 1;
            }
            if size > max_size {
                max_size = size;
            }
        }
        (self.num_clusters, singletons, non_singletons, max_size)
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
pub fn cut(tree: &Tree, method: Method) -> Result<Partition, String> {
    if tree.is_empty() {
        return Ok(Partition::new());
    }

    match method {
        Method::K(k) => simple::cut_k(tree, k),
        Method::Height(h) => simple::cut_height(tree, h),
        Method::RootDist(d) => simple::cut_root_dist(tree, d),
        Method::MaxClade(t) => clade::cut_max_clade(tree, t),
        Method::AvgClade(t) => clade::cut_avg_clade(tree, t),
        Method::MedClade(t) => clade::cut_med_clade(tree, t),
        Method::SumBranch(t) => clade::cut_sum_branch(tree, t),
        Method::Inconsistent(t, d) => inconsistent::cut_inconsistent(tree, t, d),
        Method::SingleLinkage(t) => link::cut_single_linkage(tree, t),
    }
}

// --- Helpers ---

/// Compute max distance from each node to its leaves
pub(crate) fn compute_heights(tree: &Tree, root: NodeId) -> Result<HashMap<NodeId, f64>, String> {
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
pub(crate) fn assign_clusters(
    tree: &Tree,
    cluster_roots: Vec<NodeId>,
) -> Result<Partition, String> {
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
