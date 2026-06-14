use super::Tree;
use crate::libs::phylo::node::NodeId;
use std::collections::HashMap;
use std::io::Read;

/// Read a Newick tree from a file.
///
/// # Arguments
/// * `infile` - Path to the input file (or "stdin" for stdin).
///
/// # Example
/// ```
/// // usage in CLI:
/// // let trees = pgr::libs::phylo::tree::io::from_file("path/to/tree.nwk")?;
/// ```
pub fn from_file(infile: &str) -> anyhow::Result<Vec<Tree>> {
    let mut reader = crate::reader(infile);
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
    let node = tree.get_node(node_id).unwrap();
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

/// Serialize the tree to a Graphviz DOT string.
pub fn to_dot(tree: &Tree) -> String {
    let mut s = String::from("digraph Tree {\n");
    s.push_str("    node [shape=box];\n"); // Optional styling

    if let Some(root) = tree.get_root() {
        // Use a traversal to visit all reachable nodes.
        // Preorder is good.
        if let Ok(nodes) = tree.preorder(&root) {
            for &node_id in &nodes {
                let node = tree.get_node(node_id).unwrap();
                if node.deleted {
                    continue;
                }

                // 1. Define Node
                // Use NodeID as the DOT identifier
                let label = node.name.as_deref().unwrap_or("");
                let mut label_attr = format!("label=\"{}\"", label);
                if label.is_empty() {
                    label_attr = format!("label=\"{}\"", node_id);
                }

                s.push_str(&format!("    {} [{}];\n", node_id, label_attr));

                // 2. Define Edges to children
                for &child_id in &node.children {
                    let child = tree.get_node(child_id).unwrap();
                    let mut edge_attrs = Vec::new();
                    if let Some(len) = child.length {
                        edge_attrs.push(format!("label=\"{}\"", len));
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
    }

    s.push_str("}\n");
    s
}

/// Serialize tree to SVG format.
///
/// # Arguments
/// * `tree` - The tree to serialize.
/// * `height` - Tree height for scaling branch lengths. If 0.0, uses cladogram mode.
/// * `vskip` - Vertical spacing between leaf nodes in pixels.
/// * `width` - SVG width in pixels.
pub fn to_svg(tree: &Tree, height: f64, vskip: f64, width: f64) -> String {
    let root = match tree.get_root() {
        Some(r) => r,
        None => return String::new(),
    };

    let positions = compute_svg_positions(tree, root, height, vskip, width);

    // Calculate bounding box
    let mut max_x = 0.0_f64;
    let mut max_y = 0.0_f64;
    for (_, (x, y)) in &positions {
        max_x = max_x.max(*x);
        max_y = max_y.max(*y);
    }

    // Estimate the longest leaf label width (rough: ~7px per char at 12px sans-serif)
    let max_label_width = if let Ok(nodes) = tree.preorder(&root) {
        nodes
            .iter()
            .filter_map(|&id| {
                tree.get_node(id)
                    .filter(|n| n.is_leaf() && !n.deleted)
                    .and_then(|n| n.name.as_ref())
                    .map(|name| name.replace('_', " ").len() as f64 * 7.0 + 6.0)
            })
            .fold(0.0_f64, f64::max)
    } else {
        0.0
    };

    // Add margins for labels
    let margin_left = 20.0;
    let margin_right = max_label_width.max(60.0); // at least 60px for short labels
    let margin_top = 10.0;
    let margin_bottom = 30.0; // space for scale bar

    let svg_width = max_x + margin_left + margin_right;
    let svg_height = max_y + margin_top + margin_bottom;

    let mut s = String::new();

    // SVG header
    s.push_str(&format!(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
         <svg xmlns=\"http://www.w3.org/2000/svg\" \
         width=\"{}\" height=\"{}\" viewBox=\"0 0 {} {}\">\n",
        svg_width, svg_height, svg_width, svg_height
    ));

    // Embedded styles matching template.tex
    s.push_str(
        "<style>\n\
         \tline { stroke: rgb(129,130,132); stroke-width: 1pt; stroke-linecap: round; }\n\
         \ttext { font-family: sans-serif; font-size: 12px; fill: rgb(26,25,25); }\n\
         \t.dot { fill: rgb(26,25,25); }\n\
         \t.label { font-size: 10px; }\n\
         \t.scale-text { font-size: 10px; }\n\
         </style>\n",
    );

    let ox = margin_left;
    let oy = margin_top;

    // Layer 1: Draw all edges first (so labels render on top)
    if let Ok(nodes) = tree.preorder(&root) {
        for &id in &nodes {
            let node = match tree.get_node(id) {
                Some(n) => n,
                None => continue,
            };
            if node.deleted {
                continue;
            }
            let (nx, ny) = positions[&id];

            // Horizontal branch line (from parent to this node)
            if let Some(parent_id) = node.parent {
                if let Some(&(px, _)) = positions.get(&parent_id) {
                    s.push_str(&format!(
                        "\t<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\"/>\n",
                        ox + px,
                        oy + ny,
                        ox + nx,
                        oy + ny
                    ));
                }
            }

            // Vertical connector line between children
            if node.children.len() >= 2 {
                let first_child = node.children[0];
                let last_child = *node.children.last().unwrap();
                if let (Some(&(_, first_y)), Some(&(_, last_y))) =
                    (positions.get(&first_child), positions.get(&last_child))
                {
                    s.push_str(&format!(
                        "\t<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\"/>\n",
                        ox + nx,
                        oy + first_y,
                        ox + nx,
                        oy + last_y
                    ));
                }
            }
        }

        // Layer 2: Draw dots and shapes
        for &id in &nodes {
            let node = match tree.get_node(id) {
                Some(n) => n,
                None => continue,
            };
            if node.deleted {
                continue;
            }
            let (nx, ny) = positions[&id];

            if !node.is_leaf() {
                // Only draw dot if node has a name (matching Forest behavior)
                if node.name.is_some() {
                    s.push_str(&format!(
                        "\t<circle cx=\"{}\" cy=\"{}\" r=\"2\" class=\"dot\"/>\n",
                        ox + nx,
                        oy + ny
                    ));
                }
            }
        }

        // Layer 3: Draw all text labels
        for &id in &nodes {
            let node = match tree.get_node(id) {
                Some(n) => n,
                None => continue,
            };
            if node.deleted {
                continue;
            }
            let (nx, ny) = positions[&id];

            if node.is_leaf() {
                // Leaf label: right of node
                let label = node.name.as_deref().unwrap_or("").replace('_', " ");
                if !label.is_empty() {
                    s.push_str(&format!(
                        "\t<text x=\"{}\" y=\"{}\" text-anchor=\"start\" dominant-baseline=\"middle\">{}</text>\n",
                        ox + nx + 6.0,
                        oy + ny,
                        xml_escape(&label)
                    ));
                }
            } else if let Some(name) = &node.name {
                // Internal node label: left of node
                let label = name.replace('_', " ");
                s.push_str(&format!(
                    "\t<text x=\"{}\" y=\"{}\" text-anchor=\"end\" dominant-baseline=\"middle\" class=\"label\">{}</text>\n",
                    ox + nx - 6.0,
                    oy + ny,
                    xml_escape(&label)
                ));
            }
        }
    }

    // Scale bar (phylogram mode only)
    if height > 0.0 {
        let target_scale = height / 5.0;
        let magnitude = target_scale.log10().floor();
        let base = 10.0_f64.powf(magnitude);

        let scale = [1.0, 2.0, 5.0]
            .iter()
            .map(|&x| base * x)
            .filter(|&x| x <= target_scale)
            .last()
            .unwrap_or(base);

        let scale_px = scale * (max_x / height);

        let bar_x = ox;
        let bar_y = svg_height - 20.0;

        s.push_str(&format!(
            "\t<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\"/>\n",
            bar_x,
            bar_y,
            bar_x + scale_px,
            bar_y
        ));
        s.push_str(&format!(
            "\t<text x=\"{}\" y=\"{}\" text-anchor=\"middle\" class=\"scale-text\">{}</text>\n",
            bar_x + scale_px / 2.0,
            bar_y + 14.0,
            scale
        ));
    }

    s.push_str("</svg>\n");
    s
}

/// Compute (x, y) positions for all nodes in the tree.
///
/// Layout: left-to-right (root at left, leaves at right).
/// - Y: leaves are evenly spaced vertically; internal nodes are centered over their children.
/// - X: cladogram uses tier-based alignment (all leaves at same x); phylogram uses cumulative branch length * scale.
fn compute_svg_positions(
    tree: &Tree,
    root: NodeId,
    height: f64,
    vskip: f64,
    width: f64,
) -> HashMap<NodeId, (f64, f64)> {
    let mut positions = HashMap::new();

    // Get non-deleted leaves in order
    let leaves: Vec<NodeId> = if let Ok(nodes) = tree.preorder(&root) {
        nodes
            .iter()
            .filter(|&&id| {
                tree.get_node(id)
                    .map(|n| n.is_leaf() && !n.deleted)
                    .unwrap_or(false)
            })
            .copied()
            .collect()
    } else {
        return positions;
    };

    let leaf_count = leaves.len();
    if leaf_count == 0 {
        return positions;
    }

    let y_step = vskip;

    // Assign y positions to leaves
    for (i, &leaf_id) in leaves.iter().enumerate() {
        let y = y_step * (i + 1) as f64;
        positions.insert(leaf_id, (0.0, y)); // x will be set later
    }

    // Compute cumulative branch length from root for phylogram mode
    let cum_length = if height > 0.0 {
        let mut cl = HashMap::new();
        cl.insert(root, 0.0);
        if let Ok(nodes) = tree.preorder(&root) {
            for &id in &nodes {
                let node = match tree.get_node(id) {
                    Some(n) => n,
                    None => continue,
                };
                if node.deleted {
                    continue;
                }
                let parent_len = *cl.get(&id).unwrap_or(&0.0);
                for &child_id in &node.children {
                    if let Some(child) = tree.get_node(child_id) {
                        if child.deleted {
                            continue;
                        }
                        let edge = child.length.unwrap_or(0.0);
                        cl.insert(child_id, parent_len + edge);
                    }
                }
            }
        }
        cl
    } else {
        HashMap::new()
    };

    // Compute max depth for cladogram alignment
    let max_depth = if height == 0.0 {
        leaves
            .iter()
            .filter_map(|&id| {
                let d = node_depth(tree, id);
                if d > 0 {
                    Some(d)
                } else {
                    None
                }
            })
            .max()
            .unwrap_or(1)
    } else {
        0
    };

    // Scale factor: map tree width to the requested SVG width
    let hskip = if height > 0.0 {
        width / height // pixels per unit of branch length
    } else {
        width / (max_depth as f64 + 1.0) // pixels per depth level
    };

    // Post-order traversal to compute y (center over children) and x
    if let Ok(postorder) = tree.postorder(&root) {
        for &id in &postorder {
            let node = match tree.get_node(id) {
                Some(n) => n,
                None => continue,
            };
            if node.deleted {
                continue;
            }

            if node.is_leaf() {
                // x position for leaves
                let x = if height > 0.0 {
                    let cl = cum_length.get(&id).copied().unwrap_or(0.0);
                    cl * hskip
                } else {
                    // Cladogram: all leaves aligned at max_depth * hskip
                    max_depth as f64 * hskip
                };
                if let Some(pos) = positions.get_mut(&id) {
                    pos.0 = x;
                }
            } else {
                // Internal node: y = center of children, x from depth or branch length
                let children = &node.children;
                if children.is_empty() {
                    continue;
                }

                let y_sum: f64 = children
                    .iter()
                    .filter_map(|&c| positions.get(&c).map(|p| p.1))
                    .sum();
                let y_count = children.len() as f64;
                let y = y_sum / y_count;

                let x = if height > 0.0 {
                    let cl = cum_length.get(&id).copied().unwrap_or(0.0);
                    cl * hskip
                } else {
                    // Cladogram: use subtree height (like Forest's tier system)
                    // tier = max_depth - branch_depth + 1
                    let bd = branch_depth(tree, id);
                    (max_depth - bd + 1) as f64 * hskip
                };

                positions.insert(id, (x, y));
            }
        }
    }

    positions
}

/// Escape special XML characters in text content.
fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

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
    let node = tree.get_node(id).unwrap();
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

// almost all the operations in here
fn to_forest_node_props(tree: &Tree, id: NodeId, height: f64) -> String {
    let node = tree.get_node(id).unwrap();
    let depth = node_depth(tree, id);

    let mut repr = String::new();

    let mut name = node.name.clone().map(|x| x.replace('_', " "));
    let mut color: Option<String> = None;
    let mut label: Option<String> = None;

    // internal node's name will be treated as labels and place a dot there
    if !node.is_leaf() && name.is_some() {
        label = Some(name.clone().unwrap());
        name = None;
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

    if color.is_some() {
        if name.is_some() {
            repr = format!(
                ", \\color{{{}}}{{{}}}",
                color.clone().unwrap(),
                name.clone().unwrap()
            ) + &repr;
        }
        if label.is_some() && !label.clone().unwrap().is_empty() {
            repr += &format!(
                ", label=\\color{{{}}}{{{}}}",
                color.clone().unwrap(),
                label.clone().unwrap()
            );
        }
    } else {
        if name.is_some() {
            repr = format!("{{{}}},", name.clone().unwrap()) + &repr;
        }
        if label.is_some() && !label.clone().unwrap().is_empty() {
            repr += &format!(", label={{{}}}", label.clone().unwrap());
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

// max depth of this node's children
fn branch_depth(tree: &Tree, id: NodeId) -> usize {
    let self_depth = node_depth(tree, id);
    match tree.get_subtree(&id) {
        Ok(nodes) => nodes
            .iter()
            .map(|nid| node_depth(tree, *nid))
            .max()
            .unwrap_or(self_depth),
        Err(_) => self_depth,
    }
}

fn node_depth(tree: &Tree, id: NodeId) -> usize {
    let mut depth = 0usize;
    let mut curr = id;
    loop {
        let node = match tree.get_node(curr) {
            Some(n) => n,
            None => break,
        };
        if let Some(p) = node.parent {
            depth += 1;
            curr = p;
        } else {
            break;
        }
    }
    depth
}

// relative length
fn calc_length(edge: f64, height: f64) -> i32 {
    (edge * 100.0 / height).round() as i32
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
}
