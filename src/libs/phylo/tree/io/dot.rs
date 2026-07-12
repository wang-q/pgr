//! Graphviz DOT format writer.

use super::super::Tree;

/// Escape a string for safe use inside a DOT double-quoted label.
fn escape_dot_label(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}

/// Serialize the tree to a Graphviz DOT string.
pub fn to_dot(tree: &Tree) -> String {
    let mut s = String::from("digraph Tree {\n");
    s.push_str("    node [shape=box];\n"); // Optional styling

    if let Some(root) = tree.get_root() {
        let nodes = tree.preorder(&root);
        for &node_id in &nodes {
            let Some(node) = tree.get_node(node_id) else {
                continue;
            };

            // 1. Define Node
            // Use NodeID as the DOT identifier
            let label = node.name.as_deref().unwrap_or("");
            let mut label_attr = format!("label=\"{}\"", escape_dot_label(label));
            if label.is_empty() {
                label_attr = format!("label=\"{}\"", node_id);
            }

            s.push_str(&format!("    {} [{}];\n", node_id, label_attr));

            // 2. Define Edges to children
            for &child_id in &node.children {
                let Some(child) = tree.get_node(child_id) else {
                    continue;
                };
                let mut edge_attrs = Vec::new();
                if let Some(len) = child.length {
                    if len.is_finite() && len >= 0.0 {
                        edge_attrs.push(format!("label=\"{}\"", len));
                    }
                }

                let edge_attr_str = if edge_attrs.is_empty() {
                    String::new()
                } else {
                    format!(" [{}]", edge_attrs.join(", "))
                };

                s.push_str(&format!(
                    "    {} -> {}{};\n",
                    node_id, child_id, edge_attr_str
                ));
            }
        }
    }

    s.push_str("}\n");
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_to_dot() {
        let mut tree = Tree::new();
        let n0 = tree.add_node();
        let n1 = tree.add_node();

        tree.set_root(n0);
        tree.add_child(n0, n1).unwrap();

        tree.get_node_mut(n0).unwrap().set_name("Root");
        tree.get_node_mut(n1).unwrap().set_name("A");
        tree.get_node_mut(n1).unwrap().length = Some(0.1);

        let dot = to_dot(&tree);
        assert!(dot.contains("digraph Tree {"));
        assert!(dot.contains(&format!("{} [label=\"Root\"];", n0)));
        assert!(dot.contains(&format!("{} [label=\"A\"];", n1)));
        assert!(dot.contains(&format!("{} -> {} [label=\"0.1\"];", n0, n1)));
    }

    #[test]
    fn test_to_dot_negative_length() {
        let mut tree = Tree::new();
        let n0 = tree.add_node();
        let n1 = tree.add_node();

        tree.set_root(n0);
        tree.add_child(n0, n1).unwrap();

        tree.get_node_mut(n0).unwrap().set_name("Root");
        tree.get_node_mut(n1).unwrap().set_name("A");
        tree.get_node_mut(n1).unwrap().length = Some(-0.5);

        let dot = to_dot(&tree);
        // Negative length should be treated as 0.0 (no label emitted)
        assert!(!dot.contains("label=\"-0.5\""));
        // Edge should exist but without label attribute
        assert!(dot.contains(&format!("{} -> {};", n0, n1)));
    }
}
