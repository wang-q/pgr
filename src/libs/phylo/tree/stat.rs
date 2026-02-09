use super::Tree;
use crate::libs::phylo::node::NodeId;
use std::collections::{BTreeMap, HashMap, HashSet};

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
    tree.nodes.iter().all(|n| {
        n.deleted || n.children.is_empty() || n.children.len() == 2
    })
}

/// Get bipartitions (splits) induced by the tree.
/// Returns a set of BitSets (represented as sorted vector of leaf names or IDs).
/// Here we return Set<BTreeSet<String>> for comparison.
pub fn get_partitions(tree: &Tree) -> HashSet<std::collections::BTreeSet<String>> {
    let mut partitions = HashSet::new();
    let root = match tree.get_root() {
        Some(r) => r,
        None => return partitions,
    };

    // DFS to get leaves under each node
    // Using explicit stack for postorder simulation to build sets bottom-up?
    // Or just simple recursion.
    
    fn collect_leaves(
        tree: &Tree, 
        id: NodeId, 
        partitions: &mut HashSet<std::collections::BTreeSet<String>>
    ) -> std::collections::BTreeSet<String> {
        let node = tree.get_node(id).unwrap();
        if node.children.is_empty() {
            let mut set = std::collections::BTreeSet::new();
            if let Some(name) = &node.name {
                set.insert(name.clone());
            }
            return set;
        }

        let mut current_leaves = std::collections::BTreeSet::new();
        for &child in &node.children {
            let child_leaves = collect_leaves(tree, child, partitions);
            // Add non-trivial partitions (child edge split)
            if !child_leaves.is_empty() { // && child_leaves.len() < total_leaves (handled by caller context usually)
                partitions.insert(child_leaves.clone());
                current_leaves.extend(child_leaves);
            }
        }
        current_leaves
    }

    let root_leaves = collect_leaves(tree, root, &mut partitions);
    if !root_leaves.is_empty() {
        partitions.insert(root_leaves);
    }
    
    // Partitions usually imply the split of ALL leaves into (A, B).
    // The set of leaves under a node `S` defines split `{S, All \ S}`.
    // We store `S`. When comparing, we need to handle that `S` and `All \ S` are same split.
    // For RF distance, usually we normalize by storing the smaller set, or ensure rooted trees match.
    // Here we just return the sets of leaves for each branch.
    
    // Filter out trivial partitions (size 0 or size 1 or size ALL) if desired?
    // Usually RF distance excludes trivial splits (leaf edges).
    // partitions.retain(|s| s.len() > 1 && s.len() < all_leaves.len());

    partitions
}

/// Calculate Robinson-Foulds distance to another tree.
/// (A + B - 2 * Intersection)
pub fn robinson_foulds(tree: &Tree, other: &Tree) -> usize {
    let p1 = get_partitions(tree);
    let p2 = get_partitions(other);

    let intersection = p1.intersection(&p2).count();
    p1.len() + p2.len() - 2 * intersection
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

/// Get all node names (leaf and internal).
pub fn get_names(tree: &Tree) -> Vec<String> {
    tree.nodes
        .iter()
        .filter_map(|n| {
            if !n.deleted {
                n.name.clone()
            } else {
                None
            }
        })
        .collect()
}

/// Get mapping of Name -> NodeId
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

/// Get values of a property for all nodes (e.g. "S", "C").
/// Returns Map<NodeId, Value>.
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

/// Get the node ID that has the longest incoming edge.
pub fn get_node_with_longest_edge(tree: &Tree) -> Option<NodeId> {
    tree.nodes
        .iter()
        .filter(|n| !n.deleted && n.length.is_some())
        .max_by(|a, b| a.length.unwrap_or(0.0).partial_cmp(&b.length.unwrap_or(0.0)).unwrap())
        .map(|n| n.id)
}
