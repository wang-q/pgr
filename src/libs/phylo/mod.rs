//! Phylogenetic tree library: Newick parsing, manipulation, comparison, and I/O.

pub mod cmp;
pub mod error;
pub mod node;
pub mod parser;
pub mod taxonomy;
pub mod tree;

/// Tree comparison trait (RF, WRF, KF distances and splits).
pub use cmp::TreeComparison;
/// Tree-level error type.
pub use error::TreeError;
/// Node type and its lightweight ID.
pub use node::{Node, NodeId};
/// Newick label sanitizer.
pub use parser::newick_safe;
/// Taxonomy table parser and type.
pub use taxonomy::{read_taxonomy, TaxonomyTable};
/// Arena-based phylogenetic tree.
pub use tree::Tree;
