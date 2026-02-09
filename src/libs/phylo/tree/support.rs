use super::Tree;
use crate::libs::phylo::node::NodeId;
use fixedbitset::FixedBitSet;
use std::collections::HashMap;

/// Build a map from leaf name to index (0..N-1).
/// Uses the first tree to establish the mapping.
pub fn build_leaf_map(tree: &Tree) -> Result<HashMap<String, usize>, String> {
    let mut map = HashMap::new();
    let mut index = 0;
    
    // Use traversal to ensure consistent order if needed, or just iterate nodes.
    // Iterating nodes.iter() follows creation order.
    // To be safe and deterministic, let's sort the leaf names?
    // nw_support uses the order in the first tree.
    // Let's collect names then sort them to be independent of input order in the file?
    // No, usually tools respect the order of the first tree or sort.
    // Let's sort them. It makes the indices deterministic regardless of leaf order in the input file.
    
    let mut leaf_names = Vec::new();
    for node in &tree.nodes {
        if !node.deleted && node.is_leaf() {
            if let Some(name) = &node.name {
                leaf_names.push(name.clone());
            } else {
                return Err("Leaf node missing name".to_string());
            }
        }
    }
    
    // Sorting ensures that if we run this on different files with same leaves but different order, 
    // we get compatible maps? No, map is local to the run.
    // But sorting is good practice.
    leaf_names.sort();
    
    for name in leaf_names {
        if !map.contains_key(&name) {
            map.insert(name, index);
            index += 1;
        }
    }
    
    Ok(map)
}

/// Compute bitsets for all nodes in the tree.
/// Returns a map NodeId -> FixedBitSet.
pub fn compute_all_bitsets(
    tree: &Tree, 
    leaf_map: &HashMap<String, usize>
) -> Result<HashMap<NodeId, FixedBitSet>, String> {
    let num_leaves = leaf_map.len();
    let mut node_bitsets = HashMap::new();
    
    if let Some(root) = tree.get_root() {
        // Post-order: Children processed before parents
        let traversal = tree.postorder(&root).map_err(|e| e.to_string())?;
        
        for id in traversal {
            let node = tree.get_node(id).unwrap();
            let mut bitset = FixedBitSet::with_capacity(num_leaves);
            
            if node.is_leaf() {
                if let Some(name) = &node.name {
                    if let Some(&idx) = leaf_map.get(name) {
                        bitset.set(idx, true);
                    }
                    // If leaf not in map (e.g. replicate has extra leaf), we ignore it here 
                    // or treat as missing data. nw_support assumes same leaves.
                }
            } else {
                for &child in &node.children {
                    if let Some(child_bs) = node_bitsets.get(&child) {
                        bitset.union_with(child_bs);
                    }
                }
            }
            node_bitsets.insert(id, bitset);
        }
    }
    
    Ok(node_bitsets)
}

/// Count clade frequencies from a list of replicate trees.
pub fn count_clades(
    trees: &[Tree], 
    leaf_map: &HashMap<String, usize>
) -> Result<HashMap<FixedBitSet, usize>, String> {
    let mut counts = HashMap::new();
    
    // This can be parallelized with rayon
    // But FixedBitSet is not Send/Sync? It is pure data, should be.
    // HashMap is not concurrent.
    // We can use fold/reduce.
    
    // Serial version first
    for tree in trees {
        let bitsets = compute_all_bitsets(tree, leaf_map)?;
        // Only count internal nodes? Or all nodes?
        // nw_support counts ALL nodes (leaves and internal).
        // But leaves always have count = total_reps (if present).
        // Usually we only care about support for internal nodes.
        // But nw_support implementation:
        // if (is_leaf) ... node_set_add ...
        // else ... union ... add_bipart_count ...
        // It calls add_bipart_count only for ELSE block (internal nodes).
        // So it counts internal nodes only.
        
        for (id, bs) in bitsets {
            let node = tree.get_node(id).unwrap();
            if !node.is_leaf() {
                *counts.entry(bs).or_insert(0) += 1;
            }
        }
    }
    
    Ok(counts)
}
