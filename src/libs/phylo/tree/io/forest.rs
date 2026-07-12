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
    let Some(node) = tree.get_node(id) else {
        return String::new();
    };
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
    let Some(node) = tree.get_node(id) else {
        return String::new();
    };
    let depth = node_depth(tree, id);

    let mut options = String::new();

    let mut name = node.name.clone().map(|x| x.replace('_', " "));
    let mut color: Option<String> = None;
    let mut label: Option<String> = None;

    // internal node's name will be treated as labels and place a dot there
    if !node.is_leaf() && name.is_some() {
        label = name.take();
        // dot with default color
        options += ", dot";
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
                options += &format!(", {}={{{}}}", key, v.replace('_', " "));
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
            options += &format!(", comment={{{}}}", comment);
        }
    }

    if let Some(color) = &color {
        if let Some(label) = &label {
            if !label.is_empty() {
                options += &format!(", label=\\color{{{}}}{{{}}}", color, label);
            }
        }
    } else if let Some(label) = &label {
        if !label.is_empty() {
            options += &format!(", label={{{}}}", label);
        }
    }

    let content = if let Some(color) = &color {
        if let Some(name) = &name {
            format!(
                r"{{\color{{{color}}}{{{name}}}}}",
                color = color,
                name = name
            )
        } else if node.is_leaf() {
            format!(r"{{\color{{{color}}}{{~}}}}", color = color)
        } else {
            String::new()
        }
    } else if let Some(name) = &name {
        format!("{{{}}}", name)
    } else if node.is_leaf() {
        "{~}".to_string() // non-breaking space in latex
    } else {
        String::new()
    };

    if height == 0.0 {
        let tier = if node.is_leaf() {
            0
        } else {
            branch_depth(tree, id) - depth
        };
        options += &format!(", tier={}", tier);
    } else {
        let edge = node.length.unwrap_or(0.0);
        let bl = calc_length(edge, height);
        options += &format!(", l={}mm, l sep=0", bl);

        if node.is_leaf() {
            // Add an invisible node to the rightmost to occupy spaces
            options += ", [{~},tier=0,edge={draw=none}]";
        }
    }

    if content.is_empty() {
        // Strip the leading comma separator when there is no node content.
        if options.starts_with(", ") {
            options.split_off(2)
        } else if options.starts_with(',') {
            options.split_off(1)
        } else {
            options
        }
    } else {
        content + &options
    }
}

// relative length
fn calc_length(edge: f64, height: f64) -> i32 {
    let edge = super::super::finite_length(Some(edge));
    (edge * 100.0 / height).round() as i32
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::phylo::tree::Tree;
    use std::collections::BTreeMap;

    #[test]
    fn to_forest_colored_leaf() {
        let mut tree = Tree::new();
        let root = tree.add_node();
        let leaf = tree.add_node();
        tree.set_root(root);
        tree.add_child(root, leaf).unwrap();

        if let Some(node) = tree.get_node_mut(leaf) {
            node.name = Some("Leaf_A".to_string());
            let mut props = BTreeMap::new();
            props.insert("color".to_string(), "red".to_string());
            node.properties = Some(props);
        }

        let output = to_forest(&tree, 0.0);
        assert!(
            output.contains(r"{\color{red}{Leaf A}},"),
            "expected colored leaf content, got: {}",
            output
        );
        assert!(
            !output.contains(r", \color{red}{Leaf A},"),
            "unexpected leading comma in colored leaf output: {}",
            output
        );
        assert!(
            !output.contains(",,"),
            "unexpected consecutive commas in forest output: {}",
            output
        );
    }
}
