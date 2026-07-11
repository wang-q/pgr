//! Shared tree-depth helpers for SVG and Forest layout.

use super::super::Tree;
use crate::libs::phylo::node::NodeId;

/// Depth of `id` from the root (root has depth 0).
pub(super) fn node_depth(tree: &Tree, id: NodeId) -> usize {
    let mut depth = 0usize;
    let mut curr = id;
    while let Some(node) = tree.get_node(curr) {
        match node.parent {
            Some(p) => {
                depth += 1;
                curr = p;
            }
            None => break,
        }
    }
    depth
}

/// Max depth of this node's children.
pub(super) fn branch_depth(tree: &Tree, id: NodeId) -> usize {
    let self_depth = node_depth(tree, id);
    tree.get_subtree(&id)
        .iter()
        .map(|nid| node_depth(tree, *nid))
        .max()
        .unwrap_or(self_depth)
}

/// Compute scale bar (scale value, bar length in mm) for a tree of given height.
pub fn compute_scale_bar(height: f64) -> (f64, i32) {
    let target_scale = height / 5.0;
    let magnitude = target_scale.log10().floor();
    let base = 10.0_f64.powf(magnitude);
    let scale = [1.0, 2.0, 5.0]
        .iter()
        .map(|&x| base * x)
        .rfind(|&x| x <= target_scale)
        .unwrap_or(base);
    let bar_mm = (scale * 100.0 / height).round() as i32;
    (scale, bar_mm)
}
