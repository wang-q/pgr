//! A *symmetric* scoring matrix to be used for clustering.
mod condensed;
mod named;
mod scoring;
mod transform;

pub use condensed::{get_condensed_index, CondensedMatrix};
pub use named::NamedMatrix;
pub use scoring::ScoringMatrix;
pub use transform::transform_matrix;

#[cfg(test)]
mod tests;
