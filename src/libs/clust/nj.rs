use crate::libs::pairmat::NamedMatrix;
use crate::libs::phylo::tree::Tree;
use anyhow::Result;

/// Build a tree from a distance matrix using the Neighbor-Joining algorithm.
///
/// NJ (Neighbor-Joining) is a bottom-up clustering method.
/// This implementation roots the tree at the midpoint of the last edge.
pub fn nj(matrix: &NamedMatrix) -> Result<Tree> {
    let names = matrix.get_names();
    let n = names.len();

    if n == 0 {
        return Ok(Tree::new());
    }
    if n == 1 {
        let mut tree = Tree::new();
        let root = tree.add_node();
        tree.set_root(root);
        tree.get_node_mut(root).unwrap().name = Some(names[0].clone());
        return Ok(tree);
    }

    let mut tree = Tree::new();

    // Active nodes: map from matrix index to Tree NodeId
    let mut active_nodes: Vec<usize> = Vec::with_capacity(n);

    // Initialize leaves
    for name in &names {
        let id = tree.add_node();
        tree.get_node_mut(id).unwrap().name = Some(name.to_string());
        active_nodes.push(id);
    }

    // Initialize distance matrix
    // We use a HashMap for sparse updates
    let mut dists = std::collections::HashMap::new();

    for i in 0..n {
        for j in (i + 1)..n {
            let d = matrix.get(i, j) as f64;
            let id1 = active_nodes[i];
            let id2 = active_nodes[j];
            dists.insert((id1.min(id2), id1.max(id2)), d);
        }
    }

    // NJ Loop
    while active_nodes.len() > 2 {
        let k = active_nodes.len();

        // 1. Calculate net divergence r
        let mut r = std::collections::HashMap::new();
        for &id in &active_nodes {
            let mut sum_d = 0.0;
            for &other in &active_nodes {
                if id == other {
                    continue;
                }
                let key = (id.min(other), id.max(other));
                sum_d += dists.get(&key).unwrap_or(&0.0);
            }
            r.insert(id, sum_d);
        }

        // 2. Find pair with min Q
        let mut min_q = f64::MAX;
        let mut pair = (0, 0); // Indices in active_nodes

        for i in 0..k {
            for j in (i + 1)..k {
                let id1 = active_nodes[i];
                let id2 = active_nodes[j];
                let key = (id1.min(id2), id1.max(id2));
                let d = *dists.get(&key).unwrap();
                let r1 = r[&id1];
                let r2 = r[&id2];

                let q = (k as f64 - 2.0) * d - r1 - r2;

                if q < min_q {
                    min_q = q;
                    pair = (i, j);
                }
            }
        }

        // 3. Merge nodes
        let (idx1, idx2) = pair;
        let id1 = active_nodes[idx1];
        let id2 = active_nodes[idx2];

        let d12 = *dists.get(&(id1.min(id2), id1.max(id2))).unwrap();
        let r1 = r[&id1];
        let r2 = r[&id2];

        let new_node = tree.add_node();

        // Calculate branch lengths
        let len1 = 0.5 * d12 + (r1 - r2) / (2.0 * (k as f64 - 2.0));
        let len2 = d12 - len1;

        // Add children
        tree.add_child(new_node, id1)
            .map_err(|e| anyhow::anyhow!(e))?;
        tree.add_child(new_node, id2)
            .map_err(|e| anyhow::anyhow!(e))?;

        tree.get_node_mut(id1).unwrap().length = Some(len1);
        tree.get_node_mut(id2).unwrap().length = Some(len2);

        // 4. Update distances
        let mut new_dists = Vec::new();
        for (idx, &other_id) in active_nodes.iter().enumerate() {
            if idx == idx1 || idx == idx2 {
                continue;
            }

            let d1 = *dists.get(&(id1.min(other_id), id1.max(other_id))).unwrap();
            let d2 = *dists.get(&(id2.min(other_id), id2.max(other_id))).unwrap();

            let d_new = 0.5 * (d1 + d2 - d12);
            new_dists.push((other_id, d_new));
        }

        // Update active_nodes
        if idx1 > idx2 {
            active_nodes.remove(idx1);
            active_nodes.remove(idx2);
        } else {
            active_nodes.remove(idx2);
            active_nodes.remove(idx1);
        }

        active_nodes.push(new_node);
        for (other_id, d) in new_dists {
            dists.insert((new_node.min(other_id), new_node.max(other_id)), d);
        }
    }

    // Final 2 nodes
    if active_nodes.len() == 2 {
        let id1 = active_nodes[0];
        let id2 = active_nodes[1];
        let d = *dists.get(&(id1.min(id2), id1.max(id2))).unwrap();

        // Create a root node between them
        let root = tree.add_node();
        tree.set_root(root);

        tree.add_child(root, id1).map_err(|e| anyhow::anyhow!(e))?;
        tree.add_child(root, id2).map_err(|e| anyhow::anyhow!(e))?;

        let len = d / 2.0;
        tree.get_node_mut(id1).unwrap().length = Some(len);
        tree.get_node_mut(id2).unwrap().length = Some(len);
    }

    Ok(tree)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::pairmat::NamedMatrix;
    use std::io::Write;

    #[test]
    fn test_nj_simple() {
        // Matrix:
        //   A B C D
        // A 0 7 11 14
        // B 7 0 6 9
        // C 11 6 0 7
        // D 14 9 7 0

        let content = "4
A 0 7 11 14
B 7 0 6 9
C 11 6 0 7
D 14 9 7 0
";
        let filename = "test_nj.phy";
        let mut file = std::fs::File::create(filename).unwrap();
        file.write_all(content.as_bytes()).unwrap();

        let mat = NamedMatrix::from_relaxed_phylip(filename);
        std::fs::remove_file(filename).unwrap(); // Cleanup

        let tree = nj(&mat).unwrap();
        let newick = tree.to_newick();

        assert!(newick.contains("A:"));
        assert!(newick.contains("B:"));
        assert!(newick.contains("C:"));
        assert!(newick.contains("D:"));
    }
}
