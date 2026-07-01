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

pub use dot::to_dot;
pub use forest::to_forest;
pub use newick::{from_file, to_newick, to_newick_subtree, to_newick_with_format};
pub use svg::to_svg;
pub use util::compute_scale_bar;
