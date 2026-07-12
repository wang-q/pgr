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

/// Return all node names in the tree.
pub fn get_names(tree: &Tree) -> Vec<String> {
    tree.nodes
        .iter()
        .filter(|n| !n.deleted)
        .filter_map(|n| n.name.clone())
        .collect()
}

/// Return a map from node name to node ID.
pub fn get_name_id(tree: &Tree) -> BTreeMap<String, NodeId> {
    tree.nodes
        .iter()
        .filter(|n| !n.deleted)
        .filter_map(|n| n.name.as_ref().map(|name| (name.clone(), n.id)))
        .collect()
}

/// Return a map from node ID to property value for a given key.
pub fn get_property_values(tree: &Tree, key: &str) -> BTreeMap<NodeId, String> {
    tree.nodes
        .iter()
        .filter(|n| !n.deleted)
        .filter_map(|n| {
            n.properties
                .as_ref()
                .and_then(|props| props.get(key))
                .map(|v| (n.id, v.clone()))
        })
        .collect()
}

/// Return the node ID with the longest edge length.
pub fn get_node_with_longest_edge(tree: &Tree) -> Option<NodeId> {
    tree.nodes
        .iter()
        .filter(|n| !n.deleted)
        .max_by(|a, b| {
            let len_a = a.length.unwrap_or(0.0);
            let len_b = b.length.unwrap_or(0.0);
            // Treat NaN lengths as 0.0 so they are never selected as the longest edge.
            let len_a = if len_a.is_nan() { 0.0 } else { len_a };
            let len_b = if len_b.is_nan() { 0.0 } else { len_b };
            match len_a
                .partial_cmp(&len_b)
                .unwrap_or(std::cmp::Ordering::Equal)
            {
                std::cmp::Ordering::Equal => {
                    // Tie-breaking: prefer internal nodes over leaves
                    let a_internal = !a.children.is_empty();
                    let b_internal = !b.children.is_empty();
                    if a_internal && !b_internal {
                        std::cmp::Ordering::Greater
                    } else if !a_internal && b_internal {
                        std::cmp::Ordering::Less
                    } else {
                        std::cmp::Ordering::Equal
                    }
                }
                ord => ord,
            }
        })
        .map(|n| n.id)
}

/// Check if tree is binary (all internal nodes have degree 2).
pub fn is_binary(tree: &Tree) -> bool {
    tree.nodes
        .iter()
        .all(|n| n.deleted || n.children.is_empty() || n.children.len() == 2)
}

/// Check if the tree is rooted (root node has degree 2).
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
    let root = match tree.get_root() {
        Some(r) => r,
        None => return 0.0,
    };

    let get_furthest = |start: NodeId| -> (NodeId, f64) {
        let mut max_dist = 0.0;
        let mut furthest_node = start;
        let mut visited = HashMap::new();
        let mut queue = VecDeque::new();

        visited.insert(start, 0.0);
        queue.push_back(start);

        while let Some(u) = queue.pop_front() {
            let d = *visited
                .get(&u)
                .expect("internal: BFS node was inserted before pop");
            if d > max_dist {
                max_dist = d;
                furthest_node = u;
            }

            let node = tree.get_node(u).expect("internal: BFS node exists in tree");
            let mut neighbors = node.children.clone();
            if let Some(p) = node.parent {
                neighbors.push(p);
            }

            for v in neighbors {
                if let std::collections::hash_map::Entry::Vacant(e) = visited.entry(v) {
                    let weight = if weighted {
                        let v_node = tree.get_node(v).expect("internal: neighbor exists in tree");
                        let u_node = tree.get_node(u).expect("internal: BFS node exists in tree");
                        let len = if v_node.parent == Some(u) {
                            v_node.length.unwrap_or(0.0)
                        } else {
                            u_node.length.unwrap_or(0.0)
                        };
                        // Treat NaN branch lengths as 0.0 so they do not poison
                        // the diameter calculation.
                        if len.is_nan() {
                            0.0
                        } else {
                            len
                        }
                    } else {
                        1.0
                    };
                    e.insert(d + weight);
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

/// Compute height of each node (max distance to a leaf in its subtree).
pub fn compute_node_heights(tree: &Tree) -> HashMap<NodeId, f64> {
    let mut heights = HashMap::new();

    if let Some(root) = tree.get_root() {
        let post_order = super::traversal::postorder(tree, root);
        for &id in &post_order {
            if let Some(node) = tree.get_node(id) {
                if node.children.is_empty() {
                    // Leaf
                    heights.insert(id, 0.0);
                } else {
                    let mut max_h = 0.0;
                    for &child in &node.children {
                        let child_h = *heights.get(&child).unwrap_or(&0.0);
                        let edge_len = tree.get_node(child).and_then(|n| n.length).unwrap_or(0.0);
                        let h = child_h + edge_len;
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

/// Calculate inconsistency index for each node.
pub fn calculate_inconsistency(
    tree: &Tree,
    heights: &HashMap<NodeId, f64>,
    depth: usize,
) -> HashMap<NodeId, f64> {
    super::balance::calculate_inconsistency(tree, heights, depth)
}

/// Count the number of cherries (pairs of sibling leaves) in the tree.
pub fn cherries(tree: &Tree) -> usize {
    super::balance::cherries(tree)
}

/// Compute the Sackin index (sum of leaf depths).
pub fn sackin(tree: &Tree) -> Option<f64> {
    super::balance::sackin(tree)
}

/// Compute the Colless index (sum of |left - right| subtree sizes).
pub fn colless(tree: &Tree) -> Option<f64> {
    super::balance::colless(tree)
}

/// Compute average clade distances for each internal node.
pub fn compute_avg_clade_distances(tree: &Tree) -> HashMap<NodeId, f64> {
    super::balance::compute_avg_clade_distances(tree)
}

/// Compute cumulative branch-length distance from root to each node.
pub fn compute_root_distances(tree: &Tree) -> HashMap<NodeId, f64> {
    let mut dists = HashMap::new();
    if let Some(root) = tree.get_root() {
        let mut stack = vec![(root, 0.0)];
        while let Some((node_id, d)) = stack.pop() {
            dists.insert(node_id, d);
            if let Some(node) = tree.get_node(node_id) {
                for &child in &node.children {
                    let len = tree.get_node(child).and_then(|n| n.length).unwrap_or(0.0);
                    stack.push((child, d + len));
                }
            }
        }
    }
    dists
}

/// Return (min, max, avg) root-to-leaf distances. Empty tree returns (0, 0, 0).
pub fn get_leaf_depth_stats(tree: &Tree) -> (f64, f64, f64) {
    let root_dists = compute_root_distances(tree);
    let mut depths = Vec::new();
    for (id, dist) in root_dists {
        if let Some(node) = tree.get_node(id) {
            if node.children.is_empty() {
                depths.push(dist);
            }
        }
    }
    if depths.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let min = depths.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let max = depths.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let avg = depths.iter().sum::<f64>() / depths.len() as f64;
    (min, max, avg)
}

/// Coarse classification of a tree by its branch-length decoration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TreeType {
    /// No edge has a length (all edges undecorated).
    Cladogram,
    /// All (or all-but-one, for the root) edges have a length.
    Phylogram,
    /// Mixed decorated/undecorated edges.
    Neither,
}

impl TreeType {
    /// Lowercase label used in `pgr nwk stat` output.
    pub fn as_str(self) -> &'static str {
        match self {
            TreeType::Cladogram => "cladogram",
            TreeType::Phylogram => "phylogram",
            TreeType::Neither => "neither",
        }
    }
}

/// Aggregate statistics for a single tree, as reported by `pgr nwk stat`.
#[derive(Debug, Clone)]
pub struct TreeSummary {
    pub nodes: usize,
    pub leaves: usize,
    pub dichotomies: usize,
    pub leaf_labels: usize,
    pub internal_labels: usize,
    pub edges_with_length: usize,
    pub edges_without_length: usize,
    pub cherries: usize,
    pub sackin: Option<f64>,
    pub colless: Option<f64>,
    pub is_rooted: bool,
    pub tree_type: TreeType,
}

/// Compute the full `pgr nwk stat` summary for `tree`.
pub fn tree_summary(tree: &Tree) -> TreeSummary {
    let mut nodes = 0usize;
    let mut leaves = 0usize;
    let mut dichotomies = 0usize;
    let mut leaf_labels = 0usize;
    let mut internal_labels = 0usize;
    let mut edges_with_length = 0usize;
    let mut edges_without_length = 0usize;

    if let Some(root) = tree.get_root() {
        let ids = tree.preorder(&root);
        for id in ids {
            let Some(node) = tree.get_node(id) else {
                continue;
            };
            nodes += 1;
            if node.is_leaf() {
                leaves += 1;
            }
            if node.children.len() == 2 {
                dichotomies += 1;
            }
            if node.name.is_some() {
                if node.is_leaf() {
                    leaf_labels += 1;
                } else {
                    internal_labels += 1;
                }
            }
            if node.length.is_some() {
                edges_with_length += 1;
            } else {
                edges_without_length += 1;
            }
        }
    }

    let tree_type = if edges_without_length == nodes {
        TreeType::Cladogram
    } else if edges_with_length == nodes || edges_with_length == nodes.saturating_sub(1) {
        TreeType::Phylogram
    } else {
        TreeType::Neither
    };

    TreeSummary {
        nodes,
        leaves,
        dichotomies,
        leaf_labels,
        internal_labels,
        edges_with_length,
        edges_without_length,
        cherries: cherries(tree),
        sackin: sackin(tree),
        colless: colless(tree),
        is_rooted: is_rooted(tree),
        tree_type,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::phylo::tree::Tree;

    #[test]
    fn get_node_with_longest_edge_nan() {
        // B has NaN length, C has a real length; the longest edge should be C.
        let mut tree = Tree::new();
        let root = tree.add_node();
        let a = tree.add_node();
        let b = tree.add_node();
        let c = tree.add_node();
        tree.set_root(root);
        tree.add_child(root, a).unwrap();
        tree.add_child(root, b).unwrap();
        tree.add_child(root, c).unwrap();
        tree.get_node_mut(a).unwrap().set_name("A");
        tree.get_node_mut(b).unwrap().set_name("B");
        tree.get_node_mut(b).unwrap().length = Some(f64::NAN);
        tree.get_node_mut(c).unwrap().set_name("C");
        tree.get_node_mut(c).unwrap().length = Some(5.0);

        let longest = get_node_with_longest_edge(&tree);
        let c_id = tree.get_node_by_name("C").unwrap();
        assert_eq!(longest, Some(c_id));
    }
}
