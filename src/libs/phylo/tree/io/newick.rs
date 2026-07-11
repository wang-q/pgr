//! Newick format reader and writer.

use super::super::Tree;
use crate::libs::phylo::node::NodeId;
use std::io::Read;

/// Read a Newick tree from a file.
///
/// # Arguments
/// * `infile` - Path to the input file (or "stdin" for stdin).
///
/// # Example
/// ```ignore
/// // usage in CLI:
/// // let trees = pgr::libs::phylo::tree::io::from_file("path/to/tree.nwk")?;
/// ```
pub fn from_file(infile: &str) -> anyhow::Result<Vec<Tree>> {
    let mut reader = crate::reader(infile)?;
    let mut newick = String::new();
    reader
        .read_to_string(&mut newick)
        .map_err(|e| anyhow::anyhow!("Read error: {}", e))?;
    Ok(Tree::from_newick_multi(newick.as_str())?)
}

/// Serialize tree to Newick string.
pub fn to_newick(tree: &Tree) -> String {
    to_newick_with_format(tree, "")
}

/// Serialize tree to Newick string with custom formatting options.
/// Currently supports indentation (empty for single line).
pub fn to_newick_with_format(tree: &Tree, indent: &str) -> String {
    if let Some(root) = tree.get_root() {
        let mut s = to_newick_recursive(tree, root, indent, 0);
        s.push(';');
        s
    } else {
        ";".to_string()
    }
}

/// Serialize a specific subtree to a Newick string.
pub fn to_newick_subtree(tree: &Tree, root: NodeId, indent: &str) -> String {
    let mut s = to_newick_recursive(tree, root, indent, 0);
    s.push(';');
    s
}

fn to_newick_recursive(tree: &Tree, node_id: NodeId, indent: &str, depth: usize) -> String {
    let node = tree
        .get_node(node_id)
        .expect("internal: traversal only visits existing nodes");
    let is_pretty = !indent.is_empty();

    // Calculate current indentation string
    let my_indent = if is_pretty {
        indent.repeat(depth)
    } else {
        String::new()
    };

    // Format node info: Label + Length + Comment
    let mut node_info = String::new();

    if let Some(name) = &node.name {
        node_info.push_str(&quote_label(name));
    }

    if let Some(len) = node.length {
        node_info.push_str(&format!(":{}", len));
    }

    if let Some(props) = &node.properties {
        if !props.is_empty() {
            node_info.push_str("[&&NHX");
            for (k, v) in props {
                if v.is_empty() {
                    node_info.push_str(&format!(":{}", k));
                } else {
                    node_info.push_str(&format!(":{}={}", k, v));
                }
            }
            node_info.push(']');
        }
    }

    if node.children.is_empty() {
        // Leaf: Indent + NodeInfo
        format!("{}{}", my_indent, node_info)
    } else {
        // Internal node
        let children_strs: Vec<String> = node
            .children
            .iter()
            .map(|&child| to_newick_recursive(tree, child, indent, depth + 1))
            .collect();

        if is_pretty {
            // (\n children \n)NodeInfo
            format!(
                "{}(\n{}\n{}){}",
                my_indent,
                children_strs.join(",\n"),
                my_indent,
                node_info
            )
        } else {
            format!("({}){}", children_strs.join(","), node_info)
        }
    }
}

fn quote_label(label: &str) -> String {
    let needs_quote = label.chars().any(|c| "(),:;[] \t\n".contains(c));
    if needs_quote {
        format!("'{}'", label)
    } else {
        label.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_to_newick() {
        let mut tree = Tree::new();
        //    Root
        //   /    \
        //  A:0.1  B:0.2
        let n0 = tree.add_node();
        let n1 = tree.add_node();
        let n2 = tree.add_node();

        tree.set_root(n0);
        tree.add_child(n0, n1).unwrap();
        tree.add_child(n0, n2).unwrap();

        tree.get_node_mut(n0).unwrap().set_name("Root");
        tree.get_node_mut(n1).unwrap().set_name("A");
        tree.get_node_mut(n1).unwrap().length = Some(0.1);
        tree.get_node_mut(n2).unwrap().set_name("B");
        tree.get_node_mut(n2).unwrap().length = Some(0.2);

        // Compact output
        assert_eq!(to_newick(&tree), "(A:0.1,B:0.2)Root;");

        // Pretty output
        let expected_pretty = "(\n  A:0.1,\n  B:0.2\n)Root;";
        assert_eq!(to_newick_with_format(&tree, "  "), expected_pretty);
    }

    #[test]
    fn test_to_newick_complex() {
        let mut tree = Tree::new();
        //      Root
        //     /    \
        //    I1     C
        //   /  \
        //  A    B
        let root = tree.add_node();
        let i1 = tree.add_node();
        let c = tree.add_node();
        let a = tree.add_node();
        let b = tree.add_node();

        tree.set_root(root);
        tree.get_node_mut(root).unwrap().set_name("Root");

        tree.add_child(root, i1).unwrap();
        tree.add_child(root, c).unwrap();
        tree.get_node_mut(i1).unwrap().set_name("I1");
        tree.get_node_mut(c).unwrap().set_name("C");

        tree.add_child(i1, a).unwrap();
        tree.add_child(i1, b).unwrap();
        tree.get_node_mut(a).unwrap().set_name("A");
        tree.get_node_mut(b).unwrap().set_name("B");

        // Pretty output with tab indentation
        let expected = "(\n\t(\n\t\tA,\n\t\tB\n\t)I1,\n\tC\n)Root;";
        assert_eq!(to_newick_with_format(&tree, "\t"), expected);
    }

    #[test]
    fn test_to_newick_special_chars() {
        let mut tree = Tree::new();
        let n0 = tree.add_node();
        tree.set_root(n0);
        tree.get_node_mut(n0).unwrap().set_name("Homo sapiens");

        assert_eq!(to_newick(&tree), "'Homo sapiens';");

        tree.get_node_mut(n0).unwrap().set_name("func(x)");
        assert_eq!(to_newick(&tree), "'func(x)';");
    }

    #[test]
    fn test_to_newick_properties() {
        let mut tree = Tree::new();
        let n0 = tree.add_node();
        tree.set_root(n0);
        tree.get_node_mut(n0).unwrap().set_name("A");
        tree.get_node_mut(n0).unwrap().add_property("color", "red");

        let output = to_newick(&tree);
        // Since BTreeMap order is deterministic (alphabetical keys), but we only have one key here.
        assert!(output.contains("A[&&NHX:color=red];"));
    }
}
