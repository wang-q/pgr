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

pub fn get_names(tree: &Tree) -> Vec<String> {
    tree.nodes
        .iter()
        .filter(|n| !n.deleted)
        .filter_map(|n| n.name.clone())
        .collect()
}

pub fn get_name_id(tree: &Tree) -> BTreeMap<String, NodeId> {
    tree.nodes
        .iter()
        .filter(|n| !n.deleted)
        .filter_map(|n| n.name.as_ref().map(|name| (name.clone(), n.id)))
        .collect()
}

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

pub fn get_node_with_longest_edge(tree: &Tree) -> Option<NodeId> {
    tree.nodes
        .iter()
        .filter(|n| !n.deleted)
        .max_by(|a, b| {
            let len_a = a.length.unwrap_or(0.0);
            let len_b = b.length.unwrap_or(0.0);
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
                        if v_node.parent == Some(u) {
                            v_node.length.unwrap_or(0.0)
                        } else {
                            u_node.length.unwrap_or(0.0)
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
                            let len = tree.get_node(child).and_then(|n| n.length).unwrap_or(0.0);
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
            if n == 0 {
                inconsistency.insert(id, 0.0);
            } else {
                let mean = sub_heights.iter().sum::<f64>() / n as f64;
                // Use sample variance (divisor n-1) to match SciPy/standard stats
                let divisor = if n > 1 { (n - 1) as f64 } else { 1.0 };
                let variance =
                    sub_heights.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / divisor;
                let std = variance.sqrt();
                if std == 0.0 {
                    inconsistency.insert(id, 0.0);
                } else {
                    inconsistency.insert(id, (h - mean) / std);
                }
            }
        }
    }

    inconsistency
}

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
        if let Ok(ids) = tree.preorder(&root) {
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
