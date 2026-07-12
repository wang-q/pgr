//! SVG format writer.

use super::super::Tree;
use super::util::{branch_depth, compute_scale_bar, node_depth};
use crate::libs::phylo::node::NodeId;
use std::collections::HashMap;

/// Serialize tree to SVG format.
pub fn to_svg(tree: &Tree, height: f64, vskip: f64, width: f64) -> String {
    let root = match tree.get_root() {
        Some(r) => r,
        None => return String::new(),
    };

    let positions = compute_svg_positions(tree, root, height, vskip, width);

    // Calculate bounding box
    let mut max_x = 0.0_f64;
    let mut max_y = 0.0_f64;
    for (x, y) in positions.values() {
        max_x = max_x.max(*x);
        max_y = max_y.max(*y);
    }

    // Estimate the longest leaf label width (rough: ~7px per char at 12px sans-serif)
    let nodes = tree.preorder(&root);
    let max_label_width = nodes
        .iter()
        .filter_map(|&id| {
            tree.get_node(id)
                .filter(|n| n.is_leaf() && !n.deleted)
                .and_then(|n| n.name.as_ref())
                .map(|name| name.replace('_', " ").len() as f64 * 7.0 + 6.0)
        })
        .fold(0.0_f64, f64::max);

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
    let nodes = tree.preorder(&root);
    for &id in &nodes {
        let node = match tree.get_node(id) {
            Some(n) => n,
            None => continue,
        };
        if node.deleted {
            continue;
        }
        let Some(&(nx, ny)) = positions.get(&id) else {
            continue;
        };

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
            let Some(&last_child) = node.children.last() else {
                continue;
            };
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
        let Some(&(nx, ny)) = positions.get(&id) else {
            continue;
        };

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
        let Some(&(nx, ny)) = positions.get(&id) else {
            continue;
        };

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

    // Scale bar (phylogram mode only)
    if height > 0.0 {
        let (scale, _) = compute_scale_bar(height);
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
fn compute_svg_positions(
    tree: &Tree,
    root: NodeId,
    height: f64,
    vskip: f64,
    width: f64,
) -> HashMap<NodeId, (f64, f64)> {
    let mut positions = HashMap::new();

    // Get non-deleted leaves in order
    let nodes = tree.preorder(&root);
    let leaves: Vec<NodeId> = nodes
        .iter()
        .filter(|&&id| {
            tree.get_node(id)
                .map(|n| n.is_leaf() && !n.deleted)
                .unwrap_or(false)
        })
        .copied()
        .collect();

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
        let nodes = tree.preorder(&root);
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
                    let edge = super::super::finite_length(child.length);
                    cl.insert(child_id, parent_len + edge);
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
    let postorder = tree.postorder(&root);
    for &id in &postorder {
        let node = match tree.get_node(id) {
            Some(n) => n,
            None => continue,
        };

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

            let child_ys: Vec<f64> = children
                .iter()
                .filter_map(|&c| positions.get(&c).map(|p| p.1))
                .collect();
            if child_ys.is_empty() {
                continue;
            }
            let y = child_ys.iter().sum::<f64>() / child_ys.len() as f64;

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
