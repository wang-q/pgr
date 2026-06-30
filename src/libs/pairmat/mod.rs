//! A *symmetric* scoring matrix to be used for clustering.
mod condensed;
mod named;
mod output;
mod scoring;
mod transform;

pub use condensed::{get_condensed_index, CondensedMatrix};
pub use named::NamedMatrix;
pub use output::{extract_common_upper_triangle, write_phylip_matrix, write_subset, MatrixFormat};
pub use scoring::ScoringMatrix;
pub use transform::transform_matrix;

#[cfg(test)]
mod tests;
