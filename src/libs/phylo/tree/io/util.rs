//! Shared tree-depth helpers for SVG and Forest layout.

use super::super::Tree;
use crate::libs::phylo::node::NodeId;
use std::collections::HashMap;

/// Compute depth (edges from root) for all nodes in a single BFS pass.
pub(super) fn compute_depths(tree: &Tree) -> HashMap<NodeId, usize> {
    let mut depths = HashMap::new();
    if let Some(root) = tree.get_root() {
        depths.insert(root, 0);
        for &id in &tree.levelorder(&root) {
            if let Some(node) = tree.get_node(id) {
                let depth = depths.get(&id).copied().unwrap_or(0);
                for &child in &node.children {
                    depths.insert(child, depth + 1);
                }
            }
        }
    }
    depths
}

/// Compute height (max edges to any descendant leaf) for all nodes.
/// Leaves have height 0; internal nodes have height = 1 + max(child heights).
pub(super) fn compute_heights(tree: &Tree) -> HashMap<NodeId, usize> {
    let mut heights = HashMap::new();
    if let Some(root) = tree.get_root() {
        for &id in &tree.postorder(&root) {
            if let Some(node) = tree.get_node(id) {
                let h = if node.children.is_empty() {
                    0
                } else {
                    node.children
                        .iter()
                        .filter_map(|c| heights.get(c))
                        .max()
                        .map(|h| h + 1)
                        .unwrap_or(0)
                };
                heights.insert(id, h);
            }
        }
    }
    heights
}

/// Compute scale bar (scale value, bar length in mm) for a tree of given height.
pub fn compute_scale_bar(height: f64) -> (f64, i32) {
    if height <= 0.0 || !height.is_finite() {
        return (0.0, 0);
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compute_scale_bar_non_positive() {
        assert_eq!(compute_scale_bar(0.0), (0.0, 0));
        assert_eq!(compute_scale_bar(-1.0), (0.0, 0));
    }

    #[test]
    fn compute_scale_bar_non_finite() {
        assert_eq!(compute_scale_bar(f64::NAN), (0.0, 0));
        assert_eq!(compute_scale_bar(f64::INFINITY), (0.0, 0));
        assert_eq!(compute_scale_bar(f64::NEG_INFINITY), (0.0, 0));
    }
}
