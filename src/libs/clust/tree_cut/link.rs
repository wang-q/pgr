use super::Partition;
use crate::libs::phylo::tree::Tree;
use std::collections::HashMap;

/// Cut tree using Single Linkage (cut long branches).
pub fn cut_single_linkage(tree: &Tree, threshold: f64) -> Result<Partition, String> {
    let root = tree.get_root().ok_or("Tree has no root")?;
    let mut part = Partition::new();
    let mut next_cluster_id = 0; // Starts from 0, incremented before use

    // Map NodeId -> ClusterId
    // This map is needed because when we assign a new cluster ID to a node,
    // we might not know if it's a leaf.
    // But Partition only stores Leaf assignments.

    // Stack: (node_id, cluster_id)
    // Root always starts a new cluster
    next_cluster_id += 1;
    let mut stack = vec![(root, next_cluster_id)];

    while let Some((u, cid)) = stack.pop() {
        let node = tree.get_node(u).unwrap();

        if node.children.is_empty() {
            part.assignment.insert(u, cid);
        } else {
            for &v in &node.children {
                let child_node = tree.get_node(v).unwrap();
                let len = child_node.length.unwrap_or(0.0);

                if len > threshold {
                    // Cut! v starts new cluster
                    next_cluster_id += 1;
                    stack.push((v, next_cluster_id));
                } else {
                    // v continues u's cluster
                    stack.push((v, cid));
                }
            }
        }
    }

    // Renumber clusters to be contiguous 1..K
    // Because "next_cluster_id" might have gaps if a cluster has no leaves?
    // Actually, with this logic, every created cluster ID is assigned to a node.
    // If that node is a leaf, it gets into partition.
    // If that node is internal but all its children are cut away, it has no leaves?
    // Yes, an internal node could be a cluster by itself but contain no leaves in partition map.
    // e.g. Root -> (len>T) Child. Root is cluster 1. Child is cluster 2.
    // Root has no leaves directly attached? If Root is internal node.
    // Then Cluster 1 is empty in terms of leaves.

    // So we need to normalize cluster IDs based on actual leaf assignments.
    let mut old_to_new = HashMap::new();
    let mut new_id_counter = 0;

    for val in part.assignment.values_mut() {
        if let Some(&new_id) = old_to_new.get(val) {
            *val = new_id;
        } else {
            new_id_counter += 1;
            old_to_new.insert(*val, new_id_counter);
            *val = new_id_counter;
        }
    }
    part.num_clusters = new_id_counter;

    Ok(part)
}
