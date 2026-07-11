//! LaTeX Forest format writer.

use super::super::Tree;
use super::util::{branch_depth, node_depth};
use crate::libs::phylo::node::NodeId;

/// Serialize tree to LaTeX Forest format.
///
/// # Arguments
/// * `tree` - The tree to serialize.
/// * `height` - Tree height for scaling branch lengths. If 0.0, uses cladogram mode (tier-based).
pub fn to_forest(tree: &Tree, height: f64) -> String {
    if let Some(root) = tree.get_root() {
        to_forest_recursive(tree, root, height)
    } else {
        String::new()
    }
}

fn to_forest_recursive(tree: &Tree, id: NodeId, height: f64) -> String {
    let node = tree
        .get_node(id)
        .expect("internal: traversal only visits existing nodes");
    let indent = "  ";

    let children = &node.children;
    let depth = node_depth(tree, id);

    if children.is_empty() {
        let indention = indent.repeat(depth);
        format!(
            "{}[{}]\n",
            indention,
            to_forest_node_props(tree, id, height)
        )
    } else {
        let branch_set = children
            .iter()
            .map(|&child| to_forest_recursive(tree, child, height))
            .collect::<Vec<_>>();

        let indention = indent.repeat(depth);
        format!(
            "{}[{}\n{}{}]\n",
            indention,
            to_forest_node_props(tree, id, height),
            branch_set.join(""),
            indention,
        )
    }
}

fn to_forest_node_props(tree: &Tree, id: NodeId, height: f64) -> String {
    let node = tree
        .get_node(id)
        .expect("internal: caller ensures node exists");
    let depth = node_depth(tree, id);

    let mut repr = String::new();

    let mut name = node.name.clone().map(|x| x.replace('_', " "));
    let mut color: Option<String> = None;
    let mut label: Option<String> = None;

    // internal node's name will be treated as labels and place a dot there
    if !node.is_leaf() && name.is_some() {
        label = name.take();
        // dot with default color
        repr += ", dot";
    }

    if let Some(props) = node.properties.as_ref() {
        if let Some(v) = props.get("color") {
            color = Some(v.replace('_', " "));
        }
        if let Some(v) = props.get("label") {
            label = Some(v.replace('_', " "));
        }
        for key in ["dot", "bar", "rec", "tri"] {
            if let Some(v) = props.get(key) {
                repr += &format!(", {}={{{}}}", key, v.replace('_', " "));
            }
        }
        let mut comment = String::new();
        for key in ["comment", "T", "S", "rank", "member"] {
            if let Some(v) = props.get(key) {
                if !comment.is_empty() {
                    comment += " ";
                }
                comment += &v.replace('_', " ");
            }
        }
        if !comment.is_empty() && node.is_leaf() {
            repr += &format!(", comment={{{}}}", comment);
        }
    }

    if let Some(color) = &color {
        if let Some(name) = &name {
            repr = format!(", \\color{{{}}}{{{}}}", color, name) + &repr;
        }
        if let Some(label) = &label {
            if !label.is_empty() {
                repr += &format!(", label=\\color{{{}}}{{{}}}", color, label);
            }
        }
    } else {
        if let Some(name) = &name {
            repr = format!("{{{}}},", name) + &repr;
        }
        if let Some(label) = &label {
            if !label.is_empty() {
                repr += &format!(", label={{{}}}", label);
            }
        }
    }

    if name.is_none() {
        if node.is_leaf() {
            repr = "{~},".to_owned() + &repr; // non-breaking space in latex
        } else {
            repr = ",".to_owned() + &repr;
        }
    }

    if height == 0.0 {
        let tier = if node.is_leaf() {
            0
        } else {
            branch_depth(tree, id) - depth
        };
        repr += &format!(", tier={}", tier);
    } else {
        let edge = node.length.unwrap_or(0.0);
        let bl = calc_length(edge, height);
        repr += &format!(", l={}mm, l sep=0", bl);

        if node.is_leaf() {
            // Add an invisible node to the rightmost to occupy spaces
            repr += ", [{~},tier=0,edge={draw=none}]";
        }
    }

    repr
}

// relative length
fn calc_length(edge: f64, height: f64) -> i32 {
    (edge * 100.0 / height).round() as i32
}
