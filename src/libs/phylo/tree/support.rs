use super::Tree;
use crate::libs::phylo::node::NodeId;
use fixedbitset::FixedBitSet;
use std::collections::{BTreeMap, HashMap};

/// Build a map from leaf name to index (0..N-1).
/// Uses the first tree to establish the mapping.
pub fn build_leaf_map(tree: &Tree) -> anyhow::Result<BTreeMap<String, usize>> {
    let mut map = BTreeMap::new();
    let mut index = 0;

    let mut leaf_names = Vec::new();
    for node in &tree.nodes {
        if !node.deleted && node.is_leaf() {
            if let Some(name) = &node.name {
                leaf_names.push(name.clone());
            } else {
                anyhow::bail!("Leaf node missing name");
            }
        }
    }

    leaf_names.sort();
    let len_before = leaf_names.len();
    leaf_names.dedup();
    let len_after = leaf_names.len();
    if len_after < len_before {
        log::warn!(
            "{} duplicate leaf name(s) ignored in leaf map",
            len_before - len_after
        );
    }

    for name in leaf_names {
        map.insert(name, index);
        index += 1;
    }

    Ok(map)
}

/// Compute bitsets for all nodes in the tree.
/// Returns a map NodeId -> FixedBitSet.
pub fn compute_all_bitsets(
    tree: &Tree,
    leaf_map: &BTreeMap<String, usize>,
) -> anyhow::Result<HashMap<NodeId, FixedBitSet>> {
    let num_leaves = leaf_map.len();
    let mut node_bitsets = HashMap::new();

    if let Some(root) = tree.get_root() {
        let traversal = tree.postorder(&root);

        for id in traversal {
            let Some(node) = tree.get_node(id) else {
                continue;
            };
            let mut bitset = FixedBitSet::with_capacity(num_leaves);

            if node.is_leaf() {
                if let Some(name) = &node.name {
                    if let Some(&idx) = leaf_map.get(name) {
                        bitset.set(idx, true);
                    }
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

/// Annotate internal nodes of `target` with support values from `counts`.
/// If `as_percent` is true, values are written as integer percentages of `total_reps`.
pub fn annotate_support(
    target: &mut Tree,
    leaf_map: &BTreeMap<String, usize>,
    counts: &HashMap<FixedBitSet, usize>,
    total_reps: usize,
    as_percent: bool,
) -> anyhow::Result<()> {
    let target_bitsets = compute_all_bitsets(target, leaf_map)?;
    for (id, bs) in target_bitsets {
        let node = target
            .get_node_mut(id)
            .ok_or_else(|| anyhow::anyhow!("invalid node id"))?;
        if !node.is_leaf() {
            let count = counts.get(&bs).copied().unwrap_or(0);
            let label = if as_percent {
                match (count * 100).checked_div(total_reps) {
                    Some(v) => format!("{}", v),
                    None => "0".to_string(),
                }
            } else {
                format!("{}", count)
            };
            node.name = Some(label);
        }
    }
    Ok(())
}

/// Count clade frequencies from a list of replicate trees.
pub fn count_clades(
    trees: &[Tree],
    leaf_map: &BTreeMap<String, usize>,
) -> anyhow::Result<HashMap<FixedBitSet, usize>> {
    let mut counts = HashMap::new();

    for tree in trees {
        let bitsets = compute_all_bitsets(tree, leaf_map)?;

        for (id, bs) in bitsets {
            let Some(node) = tree.get_node(id) else {
                continue;
            };
            if !node.is_leaf() {
                *counts.entry(bs).or_insert(0) += 1;
            }
        }
    }

    Ok(counts)
}
