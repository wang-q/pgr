use super::LabelMap;
use crate::libs::clust::feature::FeatureVector;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Represents a set of coordinates for items: Item -> Vector
#[derive(Debug, Clone)]
pub struct Coordinates {
    pub data: HashMap<String, Vec<f64>>,
    pub dim: usize,
}

impl Coordinates {
    /// Load coordinates from a FeatureVector file.
    /// Format: Name `tab` Val1,Val2,Val3...
    pub fn from_path<P: AsRef<Path>>(path: P) -> anyhow::Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut data = HashMap::new();
        let mut dim = 0;

        for (i, line) in reader.lines().enumerate() {
            let line = line?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }

            // Only support FeatureVector format: name \t v1,v2,v3
            let fv = FeatureVector::parse(&line)?;
            if !fv.name().is_empty() {
                let vec: Vec<f64> = fv.list().iter().map(|&v| v as f64).collect();
                if i == 0 {
                    dim = vec.len();
                } else if vec.len() != dim {
                    return Err(anyhow::anyhow!(
                        "Inconsistent dimensions at line {}: expected {}, got {}",
                        i + 1,
                        dim,
                        vec.len()
                    ));
                }
                data.insert(fv.name().clone(), vec);
            } else {
                return Err(anyhow::anyhow!(
                    "Invalid FeatureVector format at line {}: expected 'Name<tab>Val1,Val2...'",
                    i + 1
                ));
            }
        }

        if data.is_empty() {
            return Err(anyhow::anyhow!("No coordinate data found"));
        }

        Ok(Coordinates { data, dim })
    }
}

/// Calculate Davies-Bouldin Index
///
/// The score is defined as the average similarity measure of each cluster with its most similar cluster,
/// where similarity is the ratio of within-cluster distances to between-cluster distances.
///
/// The minimum score is zero, with lower values indicating better clustering.
pub fn davies_bouldin_score(partition: &LabelMap, coords: &Coordinates) -> f64 {
    // 1. Group items by cluster
    let mut clusters: HashMap<u32, Vec<&String>> = HashMap::new();
    for (item, &cluster_id) in partition {
        if coords.data.contains_key(item) {
            clusters.entry(cluster_id).or_default().push(item);
        }
    }

    let n_clusters = clusters.len();
    if n_clusters < 2 {
        // DB Index is undefined for < 2 clusters. Return 0.0 for consistency.
        return 0.0;
    }

    // 2. Calculate Centroids and Scatter (Average Intra-cluster Distance)
    struct ClusterStat {
        centroid: Vec<f64>,
        scatter: f64,
    }

    let mut stats: HashMap<u32, ClusterStat> = HashMap::new();

    for (&cluster_id, members) in &clusters {
        let n_members = members.len();
        if n_members == 0 {
            continue;
        }

        // Calculate Centroid
        let mut centroid = vec![0.0; coords.dim];
        for item in members {
            let vec = coords.data.get(*item).expect("contains_key filtered");
            for (d, val) in vec.iter().enumerate() {
                centroid[d] += val;
            }
        }
        for val in centroid.iter_mut() {
            *val /= n_members as f64;
        }

        // Calculate Scatter (Average distance to centroid)
        let mut sum_dist = 0.0;
        for item in members {
            let vec = coords.data.get(*item).expect("contains_key filtered");
            sum_dist += euclidean_dist(vec, &centroid);
        }
        let scatter = sum_dist / n_members as f64;

        stats.insert(cluster_id, ClusterStat { centroid, scatter });
    }

    // 3. Calculate DB Index
    // DB = 1/k * sum_{i} max_{j!=i} ( (s_i + s_j) / d(c_i, c_j) )
    let mut total_db = 0.0;

    for &i in stats.keys() {
        let mut max_r = 0.0; // R_ij is always >= 0
        let stat_i = stats.get(&i).expect("stat present (key from stats.keys())");

        for &j in stats.keys() {
            if i == j {
                continue;
            }
            let stat_j = stats.get(&j).expect("stat present (key from stats.keys())");

            let dist_centroids = euclidean_dist(&stat_i.centroid, &stat_j.centroid);

            // If centroids overlap perfectly, R_ij is undefined (infinity).
            // We use a large number proxy if scatter > 0, or 0.0 if scatter is also 0.
            let r_ij = if dist_centroids == 0.0 {
                if stat_i.scatter + stat_j.scatter == 0.0 {
                    0.0
                } else {
                    1e10 // Large number proxy for infinity
                }
            } else {
                (stat_i.scatter + stat_j.scatter) / dist_centroids
            };

            if r_ij > max_r {
                max_r = r_ij;
            }
        }
        total_db += max_r;
    }

    total_db / stats.len() as f64
}

/// Calculate Calinski-Harabasz Index
///
/// The score is defined as ratio between the within-cluster dispersion and
/// the between-cluster dispersion.
///
/// CH = (BGSS / (K - 1)) / (WGSS / (N - K))
///
/// Where:
/// - BGSS: Between-Group Sum of Squares (dispersion of cluster centroids from global centroid)
/// - WGSS: Within-Group Sum of Squares (dispersion of points from their cluster centroids)
///
/// Higher values indicate better clustering.
pub fn calinski_harabasz_score(partition: &LabelMap, coords: &Coordinates) -> f64 {
    // 1. Group items by cluster
    let mut clusters: HashMap<u32, Vec<&String>> = HashMap::new();
    let mut all_items = Vec::new();
    for (item, &cluster_id) in partition {
        if coords.data.contains_key(item) {
            clusters.entry(cluster_id).or_default().push(item);
            all_items.push(item);
        }
    }

    let n_clusters = clusters.len();
    let n_samples = all_items.len();

    if n_clusters < 2 || n_samples <= n_clusters {
        return 0.0;
    }

    // 2. Calculate Global Centroid
    let mut global_centroid = vec![0.0; coords.dim];
    for item in &all_items {
        let vec = coords.data.get(*item).expect("contains_key filtered");
        for (d, val) in vec.iter().enumerate() {
            global_centroid[d] += val;
        }
    }
    for val in global_centroid.iter_mut() {
        *val /= n_samples as f64;
    }

    // 3. Calculate BGSS and WGSS
    let mut bgss = 0.0;
    let mut wgss = 0.0;

    for members in clusters.values() {
        let n_members = members.len();
        if n_members == 0 {
            continue;
        }

        // Cluster Centroid
        let mut cluster_centroid = vec![0.0; coords.dim];
        for item in members {
            let vec = coords.data.get(*item).expect("contains_key filtered");
            for (d, val) in vec.iter().enumerate() {
                cluster_centroid[d] += val;
            }
        }
        for val in cluster_centroid.iter_mut() {
            *val /= n_members as f64;
        }

        // BGSS contribution: n_k * ||C_k - C||^2
        let dist_to_global = euclidean_dist(&cluster_centroid, &global_centroid);
        bgss += n_members as f64 * dist_to_global.powi(2);

        // WGSS contribution: sum ||x - C_k||^2
        for item in members {
            let vec = coords.data.get(*item).expect("contains_key filtered");
            wgss += euclidean_dist(vec, &cluster_centroid).powi(2);
        }
    }

    if wgss == 0.0 {
        if bgss == 0.0 {
            1.0 // Single point?
        } else {
            f64::INFINITY // Perfect compact clusters
        }
    } else {
        (bgss / (n_clusters as f64 - 1.0)) / (wgss / (n_samples as f64 - n_clusters as f64))
    }
}

/// Calculate PBM Index (Pakhira, Bandyopadhyay, Maulik, 2004)
///
/// PBM = ( 1/K * E_T / E_W * D_B )^2
///
/// Where:
/// - K: Number of clusters
/// - E_T: Sum of distances of all points to the global centroid (Total Scatter)
/// - E_W: Sum of within-cluster distances to centroids (Within-Group Scatter)
/// - D_B: Max distance between cluster centroids
///
/// Higher values indicate better clustering.
pub fn pbm_score(partition: &LabelMap, coords: &Coordinates) -> f64 {
    // 1. Group items by cluster
    let mut clusters: HashMap<u32, Vec<&String>> = HashMap::new();
    let mut all_items = Vec::new();
    for (item, &cluster_id) in partition {
        if coords.data.contains_key(item) {
            clusters.entry(cluster_id).or_default().push(item);
            all_items.push(item);
        }
    }

    let k = clusters.len();
    let n_samples = all_items.len();

    if k < 2 || n_samples <= k {
        return 0.0;
    }

    // 2. Calculate Global Centroid
    let mut global_centroid = vec![0.0; coords.dim];
    for item in &all_items {
        let vec = coords.data.get(*item).expect("contains_key filtered");
        for (d, val) in vec.iter().enumerate() {
            global_centroid[d] += val;
        }
    }
    for val in global_centroid.iter_mut() {
        *val /= n_samples as f64;
    }

    // 3. Calculate E_T (Total Scatter)
    let mut e_t = 0.0;
    for item in &all_items {
        let vec = coords.data.get(*item).expect("contains_key filtered");
        e_t += euclidean_dist(vec, &global_centroid);
    }

    // 4. Calculate E_W and Centroids
    let mut e_w = 0.0;
    let mut centroids: Vec<Vec<f64>> = Vec::with_capacity(k);

    for members in clusters.values() {
        let n_members = members.len();
        if n_members == 0 {
            continue;
        }

        // Cluster Centroid
        let mut cluster_centroid = vec![0.0; coords.dim];
        for item in members {
            let vec = coords.data.get(*item).expect("contains_key filtered");
            for (d, val) in vec.iter().enumerate() {
                cluster_centroid[d] += val;
            }
        }
        for val in cluster_centroid.iter_mut() {
            *val /= n_members as f64;
        }

        // E_W contribution
        for item in members {
            let vec = coords.data.get(*item).expect("contains_key filtered");
            e_w += euclidean_dist(vec, &cluster_centroid);
        }

        centroids.push(cluster_centroid);
    }

    if e_w == 0.0 {
        return if e_t > 0.0 { f64::INFINITY } else { 0.0 };
    }

    // 5. Calculate D_B (Max distance between centroids)
    let mut d_b = 0.0;
    for i in 0..centroids.len() {
        for j in i + 1..centroids.len() {
            let d = euclidean_dist(&centroids[i], &centroids[j]);
            if d > d_b {
                d_b = d;
            }
        }
    }

    // PBM = ( 1/K * E_T / E_W * D_B )^2
    let term = (1.0 / k as f64) * (e_t / e_w) * d_b;
    term * term
}

/// Calculate Ball-Hall Index
///
/// BH = 1/K * sum_{k} (1/n_k * sum_{i in C_k} ||M_i - G_k||^2)
///
/// This is the mean of the mean dispersion of each cluster.
/// Lower values indicate more compact clusters.
/// Often used with Elbow method (looking for large difference).
pub fn ball_hall_score(partition: &LabelMap, coords: &Coordinates) -> f64 {
    // 1. Group items by cluster
    let mut clusters: HashMap<u32, Vec<&String>> = HashMap::new();
    for (item, &cluster_id) in partition {
        if coords.data.contains_key(item) {
            clusters.entry(cluster_id).or_default().push(item);
        }
    }

    let k = clusters.len();
    if k == 0 {
        return 0.0;
    }

    let mut sum_mean_dispersion = 0.0;

    for members in clusters.values() {
        let n_members = members.len();
        if n_members == 0 {
            continue;
        }

        // Cluster Centroid
        let mut cluster_centroid = vec![0.0; coords.dim];
        for item in members {
            let vec = coords.data.get(*item).expect("contains_key filtered");
            for (d, val) in vec.iter().enumerate() {
                cluster_centroid[d] += val;
            }
        }
        for val in cluster_centroid.iter_mut() {
            *val /= n_members as f64;
        }

        // Sum of squared distances
        let mut sum_sq_dist = 0.0;
        for item in members {
            let vec = coords.data.get(*item).expect("contains_key filtered");
            sum_sq_dist += euclidean_dist(vec, &cluster_centroid).powi(2);
        }

        sum_mean_dispersion += sum_sq_dist / n_members as f64;
    }

    sum_mean_dispersion / k as f64
}

fn euclidean_dist(v1: &[f64], v2: &[f64]) -> f64 {
    v1.iter()
        .zip(v2.iter())
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Calculate Xie-Beni Index (Xie & Beni, 1991)
///
/// XB = WGSS / (N * min_dist_centroids^2)
///
/// Where:
/// - WGSS: Within-group sum of squares (Compactness)
/// - min_dist_centroids^2: Minimum squared Euclidean distance between cluster centers (Separation)
/// - N: Total number of points
///
/// Lower values indicate better clustering.
pub fn xie_beni_score(partition: &LabelMap, coords: &Coordinates) -> f64 {
    // 1. Group items by cluster
    let mut clusters: HashMap<u32, Vec<&String>> = HashMap::new();
    let mut all_items = Vec::new();
    for (item, &cluster_id) in partition {
        if coords.data.contains_key(item) {
            clusters.entry(cluster_id).or_default().push(item);
            all_items.push(item);
        }
    }

    let k = clusters.len();
    let n_samples = all_items.len();

    if k < 2 || n_samples <= k {
        return 0.0;
    }

    // 2. Calculate WGSS and Centroids
    let mut wgss = 0.0;
    let mut centroids: Vec<Vec<f64>> = Vec::with_capacity(k);

    for members in clusters.values() {
        let n_members = members.len();
        if n_members == 0 {
            continue;
        }

        // Cluster Centroid
        let mut cluster_centroid = vec![0.0; coords.dim];
        for item in members {
            let vec = coords.data.get(*item).expect("contains_key filtered");
            for (d, val) in vec.iter().enumerate() {
                cluster_centroid[d] += val;
            }
        }
        for val in cluster_centroid.iter_mut() {
            *val /= n_members as f64;
        }

        // WGSS contribution
        for item in members {
            let vec = coords.data.get(*item).expect("contains_key filtered");
            wgss += euclidean_dist(vec, &cluster_centroid).powi(2);
        }

        centroids.push(cluster_centroid);
    }

    // 3. Calculate Min Squared Distance between Centroids
    let mut min_sq_dist = f64::INFINITY;
    for i in 0..centroids.len() {
        for j in i + 1..centroids.len() {
            let d2 = euclidean_dist(&centroids[i], &centroids[j]).powi(2);
            if d2 < min_sq_dist {
                min_sq_dist = d2;
            }
        }
    }

    if min_sq_dist == 0.0 {
        return f64::INFINITY;
    }

    // XB = WGSS / (N * min_sq_dist)
    wgss / (n_samples as f64 * min_sq_dist)
}

/// Calculate Wemmert-Gancarski Index (Wemmert & Gancarski, 2002)
///
/// J = sum_{k=1}^K (n_k / N) * J_k
/// J_k = max(0, 1 - mean(R_k))
/// R_k(x) = ||x - G_k|| / min_{j!=k} ||x - G_j||
///
/// Measures compactness relative to separation for each point.
/// Higher values indicate better clustering (range [0, 1]).
pub fn wemmert_gancarski_score(partition: &LabelMap, coords: &Coordinates) -> f64 {
    // 1. Group items by cluster
    let mut clusters: HashMap<u32, Vec<&String>> = HashMap::new();
    let mut all_items = Vec::new();
    for (item, &cluster_id) in partition {
        if coords.data.contains_key(item) {
            clusters.entry(cluster_id).or_default().push(item);
            all_items.push(item);
        }
    }

    let k = clusters.len();
    let n_samples = all_items.len();

    if k < 2 || n_samples <= k {
        return 0.0;
    }

    // 2. Calculate Centroids
    // Map cluster ID to centroid vector
    let mut centroids: HashMap<u32, Vec<f64>> = HashMap::new();

    for (&cid, members) in &clusters {
        let n_members = members.len();
        if n_members == 0 {
            continue;
        }
        let mut cluster_centroid = vec![0.0; coords.dim];
        for item in members {
            let vec = coords.data.get(*item).expect("contains_key filtered");
            for (d, val) in vec.iter().enumerate() {
                cluster_centroid[d] += val;
            }
        }
        for val in cluster_centroid.iter_mut() {
            *val /= n_members as f64;
        }
        centroids.insert(cid, cluster_centroid);
    }

    // 3. Calculate J
    let mut weighted_j_sum = 0.0;

    for (&cid, members) in &clusters {
        let n_members = members.len();
        if n_members == 0 {
            continue;
        }

        let current_centroid = centroids.get(&cid).expect("from clusters keys");
        let mut sum_r = 0.0;

        for item in members {
            let vec = coords.data.get(*item).expect("contains_key filtered");
            let dist_intra = euclidean_dist(vec, current_centroid);

            // Find min dist to other centroids
            let mut min_dist_inter = f64::INFINITY;
            for (&other_cid, other_centroid) in &centroids {
                if cid == other_cid {
                    continue;
                }
                let d = euclidean_dist(vec, other_centroid);
                if d < min_dist_inter {
                    min_dist_inter = d;
                }
            }

            if min_dist_inter > 0.0 {
                sum_r += dist_intra / min_dist_inter;
            } else {
                // Point coincides with another centroid (bad separation).
                sum_r += 10.0; // Large penalty
            }
        }

        let mean_r = sum_r / n_members as f64;
        let j_k = (1.0 - mean_r).max(0.0);

        weighted_j_sum += (n_members as f64 / n_samples as f64) * j_k;
    }

    weighted_j_sum
}
