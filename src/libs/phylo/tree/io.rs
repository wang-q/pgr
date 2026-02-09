use super::Tree;

/// Serialize tree to Newick string.
pub fn to_newick(tree: &Tree) -> String {
    crate::libs::phylo::writer::write_newick(tree)
}

/// Serialize tree to Newick string with custom formatting options.
/// Currently supports indentation (empty for single line).
pub fn to_newick_with_format(tree: &Tree, indent: &str) -> String {
    if let Some(root) = tree.get_root() {
        crate::libs::phylo::writer::write_subtree(tree, root, indent)
    } else {
        String::new()
    }
}

/// Serialize tree to Graphviz DOT format.
pub fn to_dot(tree: &Tree) -> String {
    let mut out = String::from("graph Tree {\n");
    
    for node in &tree.nodes {
        if node.deleted { continue; }
        
        let label = node.name.clone().unwrap_or_else(|| format!("{}", node.id));
        out.push_str(&format!("  {} [label=\"{}\"];\n", node.id, label));
        
        if let Some(parent) = node.parent {
            let len = node.length.unwrap_or(1.0);
            out.push_str(&format!("  {} -- {} [len={}];\n", parent, node.id, len));
        }
    }
    
    out.push_str("}\n");
    out
}
