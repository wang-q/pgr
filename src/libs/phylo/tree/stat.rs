use super::Tree;
use crate::libs::phylo::node::NodeId;
use std::collections::{BTreeMap, HashMap};

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
    let mut leaf_counts: HashMap<NodeId, usize> = HashMap::new();
    let mut colless_sum = 0;

    // Post-order traversal ensures we process children before parents
    let nodes = tree.postorder(&root).ok()?;

    for id in nodes {
        let node = tree.get_node(id)?;

        if node.children.is_empty() {
            leaf_counts.insert(id, 1);
        } else {
            let mut sum_leaves = 0;
            let mut child_leaves = Vec::new();

            for &child in &node.children {
                let c_leaves = *leaf_counts.get(&child).unwrap_or(&0);
                sum_leaves += c_leaves;
                child_leaves.push(c_leaves);
            }
            leaf_counts.insert(id, sum_leaves);

            // Since we checked is_binary, internal nodes should have 2 children
            if node.children.len() == 2 {
                let diff = (child_leaves[0] as isize - child_leaves[1] as isize).abs();
                colless_sum += diff as usize;
            }
        }
    }

    Some(colless_sum)
}

/// Get names of all nodes in the tree (that have names).
pub fn get_names(tree: &Tree) -> Vec<String> {
    tree.nodes
        .iter()
        .filter(|n| !n.deleted)
        .filter_map(|n| n.name.clone())
        .collect()
}

/// Get a map of node name to NodeId.
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

/// Get a map of NodeId to property value for a given key.
pub fn get_property_values(tree: &Tree, key: &str) -> BTreeMap<NodeId, String> {
    let mut map = BTreeMap::new();
    for node in &tree.nodes {
        if !node.deleted {
            if let Some(props) = &node.properties {
                if let Some(val) = props.get(key) {
                    map.insert(node.id, val.clone());
                }
            }
        }
    }
    map
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
