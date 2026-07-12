//! Tree I/O: Newick, DOT, SVG, and LaTeX Forest formats.
//!
//! Each format lives in its own submodule; this module re-exports the public
//! API so existing call sites (`io::from_file`, `io::to_newick`, etc.) keep
//! working unchanged.

mod dot;
mod forest;
mod newick;
mod svg;
mod util;

/// Serialize a tree to Graphviz DOT format.
pub use dot::to_dot;
/// Serialize a tree to LaTeX Forest format.
pub use forest::to_forest;
/// Newick file reader and writers.
pub use newick::{from_file, to_newick, to_newick_subtree, to_newick_with_format};
/// Serialize a tree to SVG format.
pub use svg::to_svg;
/// Compute a scale bar for tree visualizations.
pub use util::compute_scale_bar;
