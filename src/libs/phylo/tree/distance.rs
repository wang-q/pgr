//! Distance calculations between tree nodes.
use super::Tree;
use crate::libs::phylo::node::NodeId;
use anyhow::{anyhow, Result};
use std::collections::BTreeMap;
use std::io::Write;

/// Write distance from each node (in `id_of`) to the root.
pub fn dist_root<W: Write>(
    tree: &Tree,
    id_of: &BTreeMap<String, NodeId>,
    writer: &mut W,
) -> Result<()> {
    let root = tree.get_root().ok_or_else(|| anyhow!("tree has no root"))?;
    for (k, v) in id_of.iter() {
        let dist = tree.node_distance(&root, v).map_err(anyhow::Error::msg)?;
        writer.write_fmt(format_args!("{}\t{}\n", k, format_float(dist)))?;
    }
    Ok(())
}

/// Write distance from each node to its parent.
pub fn dist_parent<W: Write>(
    tree: &Tree,
    id_of: &BTreeMap<String, NodeId>,
    writer: &mut W,
) -> Result<()> {
    for (k, v) in id_of.iter() {
        let node = tree
            .get_node(*v)
            .ok_or_else(|| anyhow!("node {} not found in tree", v))?;
        let parent = match node.parent {
            Some(p) => p,
            None => {
                writer.write_fmt(format_args!("{}\t0\n", k))?;
                continue;
            }
        };
        let dist = tree.node_distance(&parent, v).map_err(anyhow::Error::msg)?;
        writer.write_fmt(format_args!("{}\t{}\n", k, format_float(dist)))?;
    }
    Ok(())
}

/// Write pairwise distances between all nodes in `id_of`.
pub fn dist_pairwise<W: Write>(
    tree: &Tree,
    id_of: &BTreeMap<String, NodeId>,
    writer: &mut W,
) -> Result<()> {
    for (k1, v1) in id_of.iter() {
        for (k2, v2) in id_of.iter() {
            let dist = tree.node_distance(v1, v2).map_err(anyhow::Error::msg)?;
            writer.write_fmt(format_args!("{}\t{}\t{}\n", k1, k2, format_float(dist)))?;
        }
    }
    Ok(())
}

/// Write distance from each node in a pair to their Lowest Common Ancestor (LCA).
pub fn dist_lca<W: Write>(
    tree: &Tree,
    id_of: &BTreeMap<String, NodeId>,
    writer: &mut W,
) -> Result<()> {
    for (k1, v1) in id_of.iter() {
        for (k2, v2) in id_of.iter() {
            let lca = tree
                .get_common_ancestor(v1, v2)
                .map_err(anyhow::Error::msg)?;
            let dist1 = tree.node_distance(&lca, v1).map_err(anyhow::Error::msg)?;
            let dist2 = tree.node_distance(&lca, v2).map_err(anyhow::Error::msg)?;
            writer.write_fmt(format_args!(
                "{}\t{}\t{}\t{}\n",
                k1,
                k2,
                format_float(dist1),
                format_float(dist2)
            ))?;
        }
    }
    Ok(())
}

/// Write a Phylip-formatted distance matrix.
pub fn dist_phylip<W: Write>(
    tree: &Tree,
    id_of: &BTreeMap<String, NodeId>,
    writer: &mut W,
) -> Result<()> {
    let names: Vec<&String> = id_of.keys().collect();
    let n = names.len();

    // Phylip header
    writer.write_fmt(format_args!("    {}\n", n))?;

    for (i, name) in names.iter().enumerate() {
        let v1 = id_of
            .get(*name)
            .ok_or_else(|| anyhow!("node {} not found in id_of", name))?;

        // Relaxed Phylip format: name followed by space, then distances.
        writer.write_fmt(format_args!("{} ", name))?;

        for (j, other_name) in names.iter().enumerate() {
            let v2 = id_of
                .get(*other_name)
                .ok_or_else(|| anyhow!("node {} not found in id_of", other_name))?;
            let dist = if i == j {
                0.0
            } else {
                tree.node_distance(v1, v2).map_err(anyhow::Error::msg)?
            };

            writer.write_fmt(format_args!(" {:.6}", dist))?;
        }
        writer.write_all(b"\n")?;
    }
    Ok(())
}

/// Format a float by rounding to 6 decimal places and stripping trailing zeros.
fn format_float(val: f64) -> String {
    let rounded = (val * 1e6).round() / 1e6;
    format!("{}", rounded)
}
