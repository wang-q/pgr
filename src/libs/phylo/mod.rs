//! Phylogenetic tree library: Newick parsing, manipulation, comparison, and I/O.

pub mod cmp;
pub mod error;
pub mod node;
pub mod parser;
pub mod taxonomy;
pub mod tree;

pub use cmp::TreeComparison;
pub use error::TreeError;
pub use node::{Node, NodeId};
pub use parser::newick_safe;
pub use taxonomy::{read_taxonomy, TaxonomyTable};
pub use tree::Tree;
