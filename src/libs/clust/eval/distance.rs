use super::LabelMap;
use crate::libs::pairmat::NamedMatrix;
use crate::libs::phylo::tree::Tree;
use std::collections::HashMap;

/// Trait for distance matrix access
pub trait DistanceMatrix {
    fn get_distance(&self, id1: &str, id2: &str) -> f64;
}

impl DistanceMatrix for NamedMatrix {
    fn get_distance(&self, id1: &str, id2: &str) -> f64 {
        self.get_by_name(id1, id2).unwrap_or(0.0) as f64
    }
}

/// Wrapper for Tree to implement DistanceMatrix
pub struct TreeDistance {
    tree: Tree,
    name_map: HashMap<String, usize>,
}

impl TreeDistance {
    pub fn new(tree: Tree) -> Self {
        let name_map: HashMap<String, usize> = tree.get_name_id().into_iter().collect();
        Self { tree, name_map }
    }
}

impl DistanceMatrix for TreeDistance {
    fn get_distance(&self, id1: &str, id2: &str) -> f64 {
        let n1 = self.name_map.get(id1);
        let n2 = self.name_map.get(id2);
        if let (Some(&id1), Some(&id2)) = (n1, n2) {
            self.tree
                .get_distance(&id1, &id2)
                .map(|(d, _)| d)
                .unwrap_or(0.0)
        } else {
            0.0
        }
    }
}

/// Calculate Silhouette Coefficient
///
/// The Silhouette Coefficient is calculated using the mean intra-cluster distance (a)
/// and the mean nearest-cluster distance (b) for each sample.
/// The Silhouette Coefficient for a sample is (b - a) / max(a, b).
/// To clarify, b is the distance between a sample and the nearest cluster that the sample is not a part of.
/// Note that Silhouette Coefficient is only defined if number of labels is 2 <= n_labels <= n_samples - 1.
///
/// This implementation follows scikit-learn's convention:
/// - s(i) = 0 if the cluster size is 1.
pub fn silhouette_score(partition: &LabelMap, dist_mat: &dyn DistanceMatrix) -> f64 {
    // 1. Group items by cluster ID for faster access
    let mut clusters: HashMap<u32, Vec<&String>> = HashMap::new();
    for (item, &cluster_id) in partition {
        clusters.entry(cluster_id).or_default().push(item);
    }

    let n_clusters = clusters.len();
    if n_clusters < 2 || n_clusters >= partition.len() {
        return 0.0;
    }

    let mut total_s = 0.0;
    let n = partition.len();

    for (item_i, &cluster_i) in partition {
        let cluster_i_members = clusters.get(&cluster_i).unwrap();

        // Sklearn convention: s(i) = 0 if |C_i| == 1
        if cluster_i_members.len() == 1 {
            continue; // s_i is 0.0, so just skip adding
        }

        // Calculate a(i): mean distance to other items in the same cluster
        let sum_dist_a: f64 = cluster_i_members
            .iter()
            .filter(|&&item_j| item_j != item_i)
            .map(|&item_j| dist_mat.get_distance(item_i, item_j))
            .sum();
        let a_i = sum_dist_a / (cluster_i_members.len() - 1) as f64;

        // Calculate b(i): min mean distance to items in other clusters
        let mut min_mean_dist_other = f64::MAX;

        for (&cluster_j, cluster_j_members) in &clusters {
            if cluster_j == cluster_i {
                continue;
            }
            let sum_dist_b: f64 = cluster_j_members
                .iter()
                .map(|&item_j| dist_mat.get_distance(item_i, item_j))
                .sum();
            let mean_dist_b = sum_dist_b / cluster_j_members.len() as f64;
            if mean_dist_b < min_mean_dist_other {
                min_mean_dist_other = mean_dist_b;
            }
        }
        let b_i = min_mean_dist_other;

        // Calculate s(i)
        let s_i = if a_i == 0.0 && b_i == 0.0 {
            0.0
        } else {
            (b_i - a_i) / a_i.max(b_i)
        };

        total_s += s_i;
    }

    total_s / n as f64
}

/// Calculate Dunn Index
///
/// D = min(inter_cluster_dist) / max(intra_cluster_diameter)
///
/// Where:
/// - inter_cluster_dist: min distance between points in different clusters (Single Linkage)
/// - intra_cluster_diameter: max distance between points in same cluster (Complete Linkage / Diameter)
///
/// Higher values indicate better clustering.
pub fn dunn_score(partition: &LabelMap, dist_mat: &dyn DistanceMatrix) -> f64 {
    // 1. Group items by cluster ID for faster access
    let mut clusters: HashMap<u32, Vec<&String>> = HashMap::new();
    for (item, &cluster_id) in partition {
        clusters.entry(cluster_id).or_default().push(item);
    }

    let n_clusters = clusters.len();
    if n_clusters < 2 {
        return 0.0;
    }

    // 2. Calculate Max Intra-cluster Diameter (max_delta)
    let mut max_diameter = 0.0;

    for members in clusters.values() {
        let mut cluster_diameter = 0.0;
        // O(n_k^2)
        for i in 0..members.len() {
            for j in i + 1..members.len() {
                let d = dist_mat.get_distance(members[i], members[j]);
                if d > cluster_diameter {
                    cluster_diameter = d;
                }
            }
        }
        if cluster_diameter > max_diameter {
            max_diameter = cluster_diameter;
        }
    }

    // 3. Calculate Min Inter-cluster Distance (min_delta)
    let mut min_inter_dist = f64::MAX;

    let cluster_ids: Vec<u32> = clusters.keys().cloned().collect();
    for i in 0..cluster_ids.len() {
        let cid_i = cluster_ids[i];
        let members_i = clusters.get(&cid_i).unwrap();

        for cid_j in cluster_ids.iter().skip(i + 1) {
            let cid_j = *cid_j;
            let members_j = clusters.get(&cid_j).unwrap();

            // Find min dist between C_i and C_j
            for item_i in members_i {
                for item_j in members_j {
                    let d = dist_mat.get_distance(item_i, item_j);
                    if d < min_inter_dist {
                        min_inter_dist = d;
                    }
                }
            }
        }
    }

    if max_diameter == 0.0 {
        if min_inter_dist > 0.0 {
            return f64::INFINITY;
        } else {
            return 0.0;
        }
    }

    min_inter_dist / max_diameter
}

/// Calculate C-Index (Hubert & Levin, 1976)
///
/// C = (S_W - S_min) / (S_max - S_min)
///
/// Where:
/// - S_W: Sum of within-cluster distances
/// - S_min: Sum of the N_W smallest distances in the entire dataset
/// - S_max: Sum of the N_W largest distances in the entire dataset
/// - N_W: Number of within-cluster pairs
///
/// Lower values indicate better clustering (0 to 1).
/// Note: This index is computationally expensive O(N^2 log N) because it requires sorting all pairwise distances.
pub fn c_index_score(partition: &LabelMap, dist_mat: &dyn DistanceMatrix) -> f64 {
    // 1. Calculate N_W and S_W
    let mut n_w = 0;
    let mut s_w = 0.0;
    let items: Vec<&String> = partition.keys().collect();
    let n = items.len();

    if n < 2 {
        return 0.0;
    }

    let mut all_distances: Vec<f64> = Vec::with_capacity(n * (n - 1) / 2);

    for i in 0..n {
        for j in i + 1..n {
            let d = dist_mat.get_distance(items[i], items[j]);
            all_distances.push(d);

            if partition[items[i]] == partition[items[j]] {
                s_w += d;
                n_w += 1;
            }
        }
    }

    if n_w == 0 {
        // C-index is not defined for singletons (no within-cluster pairs).
        return 0.0;
    }

    if n_w == all_distances.len() {
        // Single cluster containing all points.
        return 0.0;
    }

    // 2. Sort all distances to find S_min and S_max
    // Use float sorting (handle NaNs as max or ignore)
    all_distances.sort_by(|a, b| a.total_cmp(b));

    let s_min: f64 = all_distances.iter().take(n_w).sum();
    let s_max: f64 = all_distances.iter().rev().take(n_w).sum();

    if (s_max - s_min).abs() < 1e-9 {
        return 0.0;
    }

    (s_w - s_min) / (s_max - s_min)
}

/// Calculate Hubert's Gamma
///
/// Correlation between the distance matrix and the binary matrix of cluster membership.
/// Here we define the binary matrix Y as:
/// Y_ij = 0 if i and j are in the same cluster
/// Y_ij = 1 if i and j are in different clusters
///
/// High Gamma indicates that distances between different clusters are generally larger
/// than distances within the same cluster.
/// Range: [-1, 1]. Higher is better.
pub fn gamma_score(partition: &LabelMap, dist_mat: &dyn DistanceMatrix) -> f64 {
    let items: Vec<&String> = partition.keys().collect();
    let n = items.len();
    if n < 2 {
        return 0.0;
    }

    let n_t = (n * (n - 1) / 2) as f64;

    let mut sum_x = 0.0; // Sum of distances
    let mut sum_x_sq = 0.0; // Sum of squared distances
    let mut sum_y = 0.0; // Sum of binary values (count of diff pairs)
    let mut sum_xy = 0.0; // Sum of distance * binary (Sum of between-cluster distances)

    for i in 0..n {
        for j in i + 1..n {
            let d = dist_mat.get_distance(items[i], items[j]);
            let is_diff = if partition[items[i]] != partition[items[j]] {
                1.0
            } else {
                0.0
            };

            sum_x += d;
            sum_x_sq += d * d;
            sum_y += is_diff;
            if is_diff > 0.0 {
                sum_xy += d;
            }
        }
    }

    // Pearson Correlation
    // numerator = N * sum(xy) - sum(x) * sum(y)
    // denominator = sqrt(N * sum(x^2) - sum(x)^2) * sqrt(N * sum(y^2) - sum(y)^2)

    let numerator = n_t * sum_xy - sum_x * sum_y;
    let var_x = n_t * sum_x_sq - sum_x * sum_x;
    let var_y = n_t * sum_y - sum_y * sum_y;

    if var_x <= 0.0 || var_y <= 0.0 {
        return 0.0;
    }

    numerator / (var_x.sqrt() * var_y.sqrt())
}

/// Calculate Kendall's Tau (for clustering)
///
/// Measures the ordinal association between the distance matrix and the cluster membership.
/// Like Gamma, we assume Y=0 (same), Y=1 (diff).
/// Range: [-1, 1]. Higher is better.
pub fn tau_score(partition: &LabelMap, dist_mat: &dyn DistanceMatrix) -> f64 {
    let items: Vec<&String> = partition.keys().collect();
    let n = items.len();
    if n < 2 {
        return 0.0;
    }

    // 1. Collect all pairs (distance, is_diff)
    struct Pair {
        dist: f64,
        is_diff: bool,
    }

    let mut pairs = Vec::with_capacity(n * (n - 1) / 2);
    for i in 0..n {
        for j in i + 1..n {
            let d = dist_mat.get_distance(items[i], items[j]);
            let is_diff = partition[items[i]] != partition[items[j]];
            pairs.push(Pair { dist: d, is_diff });
        }
    }

    // 2. Sort by distance
    pairs.sort_by(|a, b| a.dist.total_cmp(&b.dist));

    // 3. Calculate Concordant (S+) and Discordant (S-)
    let mut s_plus = 0.0;
    let mut s_minus = 0.0;

    let mut cum_same = 0;
    let mut cum_diff = 0;

    let mut i = 0;
    while i < pairs.len() {
        // Find block of identical distances
        let mut j = i;
        while j < pairs.len() && (pairs[j].dist - pairs[i].dist).abs() < 1e-10 {
            j += 1;
        }

        // Process block i..j
        let mut block_same = 0;
        let mut block_diff = 0;

        for p in pairs.iter().take(j).skip(i) {
            if p.is_diff {
                // This pair is DIFF.
                // Compared to previous SAME pairs (which have strictly smaller distance):
                // Concordant: Previous Same < Current Diff
                s_plus += cum_same as f64;
                block_diff += 1;
            } else {
                // This pair is SAME.
                // Compared to previous DIFF pairs (which have strictly smaller distance):
                // Discordant: Previous Diff < Current Same
                s_minus += cum_diff as f64;
                block_same += 1;
            }
        }

        cum_same += block_same;
        cum_diff += block_diff;

        i = j;
    }

    let n_pairs = pairs.len() as f64;
    let n0 = n_pairs * (n_pairs - 1.0) / 2.0;

    let n_w = cum_same as f64;
    let n_b = cum_diff as f64;
    let n2 = n_w * (n_w - 1.0) / 2.0 + n_b * (n_b - 1.0) / 2.0;

    // Ties in X: Blocks of equal distance
    let mut n1 = 0.0;
    let mut i = 0;
    while i < pairs.len() {
        let mut j = i;
        while j < pairs.len() && (pairs[j].dist - pairs[i].dist).abs() < 1e-10 {
            j += 1;
        }
        let block_len = (j - i) as f64;
        n1 += block_len * (block_len - 1.0) / 2.0;
        i = j;
    }

    let denom = ((n0 - n1) * (n0 - n2)).sqrt();

    if denom == 0.0 {
        0.0
    } else {
        (s_plus - s_minus) / denom
    }
}
