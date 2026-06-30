mod coordinates;
mod distance;
pub mod format;
mod pairwise;
mod partition;

pub use coordinates::{
    ball_hall_score, calinski_harabasz_score, davies_bouldin_score, pbm_score,
    wemmert_gancarski_score, xie_beni_score, Coordinates,
};
pub use distance::{
    c_index_score, dunn_score, gamma_score, silhouette_score, tau_score, DistanceMatrix,
    TreeDistance,
};
pub use pairwise::{evaluate, Metrics};
pub use partition::{load_batch_partitions, load_partition, remove_singletons, PartitionFormat};

use std::collections::HashMap;

/// Represents a clustering partition: Item -> ClusterID
pub type LabelMap = HashMap<String, u32>;

#[cfg(test)]
mod tests;
