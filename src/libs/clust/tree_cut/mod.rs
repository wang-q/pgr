use crate::libs::phylo::node::NodeId;
use crate::libs::phylo::tree::Tree;
use std::collections::HashMap;

pub mod clade;
pub mod dynamic;
pub mod hybrid;
pub mod inconsistent;
pub mod link;
pub mod method;
pub mod partition;
pub mod simple;

pub use method::{build_method, Method, METHOD_NAMES};
pub use partition::{
    find_representative, format_clusters, format_scan_rows, partition_to_clusters, Cluster,
    Partition, RepMode,
};

use dynamic::{cutree_dynamic_tree, DynamicTreeOptions};
use hybrid::{cutree_hybrid, HybridOptions};

/// Cut-method dispatch value (clap-free). The caller constructs this from
/// `ArgMatches` and `dispatch_cut` executes the corresponding algorithm.
pub enum CutDispatch {
    /// Dynamic Tree Cut.
    DynamicTree(DynamicTreeOptions),
    /// Dynamic Hybrid Cut (with pre-loaded distance matrix).
    DynamicHybrid(HybridOptions),
    /// One of the standard `METHOD_NAMES` methods.
    Standard {
        name: &'static str,
        val: f64,
        deep: usize,
        leaf_depths: Option<(f64, f64, f64)>,
    },
}

/// Execute the cut specified by `dispatch` on `tree`. Returns the resulting
/// partition and the method name (for labeling).
pub fn dispatch_cut(
    tree: &Tree,
    dispatch: CutDispatch,
) -> anyhow::Result<(Partition, &'static str)> {
    match dispatch {
        CutDispatch::DynamicTree(opts) => {
            let p = cutree_dynamic_tree(tree, opts).map_err(|e| anyhow::anyhow!(e))?;
            Ok((p, "dynamic-tree"))
        }
        CutDispatch::DynamicHybrid(opts) => {
            let p = cutree_hybrid(tree, opts).map_err(|e| anyhow::anyhow!(e))?;
            Ok((p, "dynamic-hybrid"))
        }
        CutDispatch::Standard {
            name,
            val,
            deep,
            leaf_depths,
        } => {
            let method =
                build_method(name, val, deep, leaf_depths).map_err(|e| anyhow::anyhow!(e))?;
            let p = cut(tree, method).map_err(|e| anyhow::anyhow!(e))?;
            Ok((p, name))
        }
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
        Method::MaxClade(t) => clade::cut_max_clade(tree, t).map_err(|e| e.to_string()),
        Method::AvgClade(t) => clade::cut_avg_clade(tree, t).map_err(|e| e.to_string()),
        Method::MedClade(t) => clade::cut_med_clade(tree, t).map_err(|e| e.to_string()),
        Method::SumBranch(t) => clade::cut_sum_branch(tree, t).map_err(|e| e.to_string()),
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

/// Mask low-support internal nodes by setting branch length to infinity (TreeCluster semantics).
pub fn apply_support_filter(tree: &mut Tree, threshold: f64) {
    let len = tree.len();
    for i in 0..len {
        let should_mask = {
            if let Some(node) = tree.get_node(i) {
                if !node.children.is_empty() {
                    let support = node
                        .name
                        .as_ref()
                        .and_then(|n| n.parse::<f64>().ok())
                        .unwrap_or(100.0);
                    support < threshold
                } else {
                    false
                }
            } else {
                false
            }
        };

        if should_mask {
            if let Some(node) = tree.get_node_mut(i) {
                node.length = Some(f64::INFINITY);
            }
        }
    }
}
