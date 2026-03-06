use crate::libs::fmt::feature::FeatureVector;
use crate::libs::pairmat::NamedMatrix;
use crate::libs::phylo::tree::Tree;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Represents a clustering partition: Item -> ClusterID
pub type Partition = HashMap<String, u32>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PartitionFormat {
    Cluster,
    Pair,
    Long,
}

impl std::str::FromStr for PartitionFormat {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "cluster" => Ok(PartitionFormat::Cluster),
            "pair" => Ok(PartitionFormat::Pair),
            "long" => Ok(PartitionFormat::Long),
            _ => Err(format!("Unknown format: {}", s)),
        }
    }
}

/// Load a partition from a file.
/// Supports two formats:
/// 1. Cluster-based: Each line is a cluster, items separated by whitespace.
///    The first item is treated as the cluster representative/ID.
/// 2. Pair-based: Two columns.
///    - If 2 columns: ClusterID <tab> Item
///    - If > 2 columns: Treated as Cluster-based.
/// 3. Long-based: Treated as Batch Partition (returns empty map here, use load_batch_partitions).
pub fn load_partition<P: AsRef<Path>>(
    path: P,
    format: PartitionFormat,
) -> anyhow::Result<Partition> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if !line.trim().is_empty() && !line.starts_with('#') {
            lines.push(line);
        }
    }

    if lines.is_empty() {
        return Ok(HashMap::new());
    }

    match format {
        PartitionFormat::Cluster => parse_cluster_format(&lines),
        PartitionFormat::Pair => parse_pair_format(&lines),
        PartitionFormat::Long => Err(anyhow::anyhow!(
            "Long format is for batch processing. Use load_batch_partitions instead."
        )),
    }
}

fn parse_pair_format(lines: &[String]) -> anyhow::Result<Partition> {
    let mut partition = HashMap::new();
    // We need to map string labels to u32 IDs
    let mut label_map: HashMap<String, u32> = HashMap::new();
    let mut next_id = 0;

    for line in lines {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 {
            continue;
        }
        let label_str = parts[0];
        let item = parts[1];

        let label_id = *label_map.entry(label_str.to_string()).or_insert_with(|| {
            next_id += 1;
            next_id
        });

        partition.insert(item.to_string(), label_id);
    }
    Ok(partition)
}

fn parse_cluster_format(lines: &[String]) -> anyhow::Result<Partition> {
    let mut partition = HashMap::new();
    let mut cluster_id = 0;

    for line in lines {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.is_empty() {
            continue;
        }
        cluster_id += 1;
        for item in parts {
            partition.insert(item.to_string(), cluster_id);
        }
    }
    Ok(partition)
}

/// Load batch partitions from a file in Long format.
/// Format: GroupID <tab> ClusterID <tab> SampleID
/// Returns a list of (GroupID, Partition).
pub fn load_batch_partitions<P: AsRef<Path>>(path: P) -> anyhow::Result<Vec<(String, Partition)>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut groups: Vec<String> = Vec::new();
    let mut group_indices: HashMap<String, usize> = HashMap::new();
    let mut partitions: Vec<Partition> = Vec::new();
    
    // Per-group label mapping to handle non-numeric cluster IDs consistently
    let mut group_label_maps: Vec<HashMap<String, u32>> = Vec::new();
    let mut group_next_ids: Vec<u32> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty()
            || line.starts_with('#')
            || line.starts_with("Threshold")
            || line.starts_with("Group")
        {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue; // Skip invalid lines
        }

        let group_id = parts[0].to_string();
        let cluster_id_str = parts[1];
        let sample_id = parts[2].to_string();

        let idx = if let Some(&idx) = group_indices.get(&group_id) {
            idx
        } else {
            let idx = groups.len();
            groups.push(group_id.clone());
            group_indices.insert(group_id, idx);
            partitions.push(HashMap::new());
            group_label_maps.push(HashMap::new());
            group_next_ids.push(0);
            idx
        };

        // Map cluster ID string to unique u32 for this group
        let cluster_id = *group_label_maps[idx]
            .entry(cluster_id_str.to_string())
            .or_insert_with(|| {
                group_next_ids[idx] += 1;
                group_next_ids[idx]
            });

        partitions[idx].insert(sample_id, cluster_id);
    }

    let result = groups.into_iter().zip(partitions.into_iter()).collect();
    Ok(result)
}

#[derive(Debug, Default)]
pub struct Metrics {
    /// Adjusted Rand Index
    pub ari: f64,
    /// Adjusted Mutual Information
    pub ami: f64,
    /// Homogeneity (each cluster contains only members of a single class)
    pub homogeneity: f64,
    /// Completeness (all members of a given class are assigned to the same cluster)
    pub completeness: f64,
    /// V-Measure (harmonic mean of homogeneity and completeness)
    pub v_measure: f64,
    /// Fowlkes-Mallows Index (geometric mean of precision and recall)
    pub fmi: f64,
    /// Normalized Mutual Information
    pub nmi: f64,
    /// Mutual Information
    pub mi: f64,
    /// Rand Index
    pub ri: f64,
    /// Jaccard Index
    pub jaccard: f64,
    /// Precision (Pair-wise)
    pub precision: f64,
    /// Recall (Pair-wise)
    pub recall: f64,
}

pub fn evaluate(p1: &Partition, p2: &Partition) -> Metrics {
    // 1. Find intersection of keys
    let keys1: HashSet<_> = p1.keys().collect();
    let keys2: HashSet<_> = p2.keys().collect();
    let common_keys: Vec<&String> = keys1.intersection(&keys2).cloned().collect();

    let n = common_keys.len();
    if n == 0 {
        return Metrics::default();
    }

    // 2. Map labels to 0..A and 0..B
    // We need contiguous integers for the contingency table
    let (labels1, n_clusters_1) = normalize_labels(p1, &common_keys);
    let (labels2, n_clusters_2) = normalize_labels(p2, &common_keys);

    // 3. Build Contingency Table
    // Map (u32, u32) -> count
    let mut table: HashMap<(u32, u32), usize> = HashMap::new();
    let mut a_counts: Vec<usize> = vec![0; n_clusters_1];
    let mut b_counts: Vec<usize> = vec![0; n_clusters_2];

    for i in 0..n {
        let u = labels1[i];
        let v = labels2[i];
        *table.entry((u, v)).or_insert(0) += 1;
        a_counts[u as usize] += 1;
        b_counts[v as usize] += 1;
    }

    // 4. Calculate Concordance Matrix
    let concordance = ConcordanceMatrix::calculate(&table, &a_counts, &b_counts, n);

    // 5. Calculate Metrics
    let ari = calculate_ari_from_concordance(&concordance);
    let (homogeneity, completeness, v_measure, mi, nmi) =
        calculate_v_measure_and_mi(&table, &a_counts, &b_counts, n);
    let ami = calculate_ami(&table, &a_counts, &b_counts, n);

    // Set matching metrics
    let ri = concordance.rand_index();
    let jaccard = concordance.jaccard_index();
    let precision = concordance.precision();
    let recall = concordance.recall();
    let fmi = concordance.fowlkes_mallows_index();

    Metrics {
        ari,
        ami,
        homogeneity,
        completeness,
        v_measure,
        fmi,
        nmi,
        mi,
        ri,
        jaccard,
        precision,
        recall,
    }
}

pub struct ConcordanceMatrix {
    pub tp: f64,  // N_yy (Same cluster in both)
    pub fn_: f64, // N_yn (Same in P1, Different in P2)
    pub fp: f64,  // N_ny (Different in P1, Same in P2)
    pub tn: f64,  // N_nn (Different in both)
    pub total_pairs: f64,
}

impl ConcordanceMatrix {
    pub fn calculate(
        table: &HashMap<(u32, u32), usize>,
        a_counts: &[usize],
        b_counts: &[usize],
        n: usize,
    ) -> Self {
        fn binom2(x: usize) -> f64 {
            if x < 2 {
                0.0
            } else {
                (x as f64 * (x as f64 - 1.0)) / 2.0
            }
        }

        let tp: f64 = table.values().map(|&count| binom2(count)).sum();
        let sum_a_2: f64 = a_counts.iter().map(|&count| binom2(count)).sum();
        let sum_b_2: f64 = b_counts.iter().map(|&count| binom2(count)).sum();
        let total_pairs = binom2(n);

        let fp = sum_a_2 - tp;
        let fn_ = sum_b_2 - tp;
        let tn = total_pairs - (tp + fp + fn_);

        ConcordanceMatrix {
            tp,
            fn_,
            fp,
            tn,
            total_pairs,
        }
    }

    pub fn rand_index(&self) -> f64 {
        if self.total_pairs == 0.0 {
            0.0
        } else {
            (self.tp + self.tn) / self.total_pairs
        }
    }

    pub fn jaccard_index(&self) -> f64 {
        let denom = self.tp + self.fp + self.fn_;
        if denom == 0.0 {
            0.0
        } else {
            self.tp / denom
        }
    }

    pub fn precision(&self) -> f64 {
        let denom = self.tp + self.fp;
        if denom == 0.0 {
            0.0
        } else {
            self.tp / denom
        }
    }

    pub fn recall(&self) -> f64 {
        let denom = self.tp + self.fn_;
        if denom == 0.0 {
            0.0
        } else {
            self.tp / denom
        }
    }

    pub fn fowlkes_mallows_index(&self) -> f64 {
        let prec = self.precision();
        let rec = self.recall();
        if prec * rec == 0.0 {
            0.0
        } else {
            (prec * rec).sqrt()
        }
    }
}

fn calculate_ari_from_concordance(c: &ConcordanceMatrix) -> f64 {
    // ARI = (RI - Expected_RI) / (Max_RI - Expected_RI)
    // Numerator is also: TP - Expected_TP
    // Expected_TP = (sum_a_2 * sum_b_2) / total_pairs
    //             = (TP+FP) * (TP+FN) / total_pairs

    if c.total_pairs == 0.0 {
        return 1.0; // Trivial case (0 or 1 item)
    }

    if c.fp == 0.0 && c.fn_ == 0.0 {
        return 1.0; // Perfect match
    }

    let index = c.tp;
    let sum_a_2 = c.tp + c.fp;
    let sum_b_2 = c.tp + c.fn_;

    let expected_index = (sum_a_2 * sum_b_2) / c.total_pairs;
    let max_index = (sum_a_2 + sum_b_2) / 2.0;

    if max_index - expected_index == 0.0 {
        return 0.0;
    }

    (index - expected_index) / (max_index - expected_index)
}

fn normalize_labels(p: &Partition, keys: &[&String]) -> (Vec<u32>, usize) {
    let mut label_map = HashMap::new();
    let mut new_labels = Vec::with_capacity(keys.len());
    let mut next_id = 0;

    for key in keys {
        let original_label = p.get(*key).unwrap();
        let id = *label_map.entry(original_label).or_insert_with(|| {
            let id = next_id;
            next_id += 1;
            id
        });
        new_labels.push(id);
    }
    (new_labels, next_id as usize)
}

fn calculate_v_measure_and_mi(
    table: &HashMap<(u32, u32), usize>,
    a_counts: &[usize],
    b_counts: &[usize],
    n: usize,
) -> (f64, f64, f64, f64, f64) {
    if n == 0 {
        return (0.0, 0.0, 0.0, 0.0, 0.0);
    }
    let n_f = n as f64;

    // Entropy H(U) = - sum (a_i/N) log (a_i/N)
    let h_a = entropy(a_counts, n_f);
    let h_b = entropy(b_counts, n_f);

    // Mutual Information MI(U, V)
    // MI = sum_ij (n_ij/N) log ( (n_ij * N) / (a_i * b_j) )
    let mut mi = 0.0;
    for (&(u, v), &nij) in table {
        if nij == 0 {
            continue;
        }
        let nij_f = nij as f64;
        let ai = a_counts[u as usize] as f64;
        let bj = b_counts[v as usize] as f64;
        let term = (nij_f / n_f) * ((nij_f * n_f) / (ai * bj)).ln();
        mi += term;
    }

    // Normalized Mutual Information (NMI)
    // NMI = MI / sqrt(H(U) * H(V))
    let nmi = if h_a * h_b == 0.0 {
        0.0
    } else {
        mi / (h_a * h_b).sqrt()
    };

    // Homogeneity = 1 - H(U|V) / H(U) = MI(U,V) / H(U)
    // Completeness = 1 - H(V|U) / H(V) = MI(U,V) / H(V)

    let homogeneity = if h_a == 0.0 { 1.0 } else { mi / h_a };
    let completeness = if h_b == 0.0 { 1.0 } else { mi / h_b };

    let v_measure = if homogeneity + completeness == 0.0 {
        0.0
    } else {
        2.0 * homogeneity * completeness / (homogeneity + completeness)
    };

    (homogeneity, completeness, v_measure, mi, nmi)
}

fn entropy(counts: &[usize], n: f64) -> f64 {
    let mut h = 0.0;
    for &count in counts {
        if count == 0 {
            continue;
        }
        let p = count as f64 / n;
        h -= p * p.ln();
    }
    h
}

/// Calculate Adjusted Mutual Information (AMI)
/// AMI = (MI - E[MI]) / (mean(H(U), H(V)) - E[MI])
fn calculate_ami(
    table: &HashMap<(u32, u32), usize>,
    a_counts: &[usize],
    b_counts: &[usize],
    n: usize,
) -> f64 {
    // AMI = (MI - E[MI]) / (mean(H(U), H(V)) - E[MI])

    let n_f = n as f64;
    let h_a = entropy(a_counts, n_f);
    let h_b = entropy(b_counts, n_f);

    let mut mi = 0.0;
    for (&(u, v), &nij) in table {
        if nij == 0 {
            continue;
        }
        let nij_f = nij as f64;
        let ai = a_counts[u as usize] as f64;
        let bj = b_counts[v as usize] as f64;
        let term = (nij_f / n_f) * ((nij_f * n_f) / (ai * bj)).ln();
        mi += term;
    }

    let expected_mi = expected_mutual_info(a_counts, b_counts, n);

    let mean_h = (h_a + h_b) / 2.0;

    if h_a == 0.0 && h_b == 0.0 {
        return 1.0;
    }

    if mean_h - expected_mi == 0.0 {
        return 0.0;
    }

    (mi - expected_mi) / (mean_h - expected_mi)
}

/// Calculate Expected Mutual Information (EMI)
/// Based on Vinh et al. (2010)
/// E[MI] = sum_ij ...
/// Using the optimized summation formula.
fn expected_mutual_info(a_counts: &[usize], b_counts: &[usize], n: usize) -> f64 {
    let n_f = n as f64;
    let mut emi = 0.0;

    for &ai in a_counts {
        for &bj in b_counts {
            if ai == 0 || bj == 0 {
                continue;
            }

            // Range of nij: max(1, ai+bj-N) .. min(ai, bj)
            let start = if ai + bj > n { ai + bj - n } else { 1 };
            let end = std::cmp::min(ai, bj);

            if start > end {
                continue;
            }

            for nij in start..=end {
                // Probability of this contingency table entry under hypergeometric model
                // P(nij) = C(ai, nij) * C(N-ai, bj-nij) / C(N, bj)
                // This is probability of overlap size `nij` given margins `ai` and `bj`.

                // We sum: (nij / N) * log( (N*nij) / (ai*bj) ) * P(nij)

                let p_nij = log_hypergeometric_prob(nij, ai, bj, n).exp();
                if p_nij == 0.0 {
                    continue;
                }

                let nij_f = nij as f64;
                let ai_f = ai as f64;
                let bj_f = bj as f64;

                let term = (nij_f / n_f) * ((n_f * nij_f) / (ai_f * bj_f)).ln() * p_nij;
                emi += term;
            }
        }
    }
    emi
}

fn log_combination(n: usize, k: usize) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    if k == 0 || k == n {
        return 0.0;
    }
    lgamma(n as f64 + 1.0) - lgamma(k as f64 + 1.0) - lgamma((n - k) as f64 + 1.0)
}

fn log_hypergeometric_prob(nij: usize, ai: usize, bj: usize, n: usize) -> f64 {
    // P = C(ai, nij) * C(N-ai, bj-nij) / C(N, bj)
    // logP = logC(ai, nij) + logC(N-ai, bj-nij) - logC(N, bj)

    log_combination(ai, nij) + log_combination(n - ai, bj - nij) - log_combination(n, bj)
}

/// Lanczos approximation for log(gamma(z))
fn lgamma(x: f64) -> f64 {
    let p = [
        676.5203681218851,
        -1259.1392167224028,
        771.323_428_777_653,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.984369578019572e-6,
        1.5056327351493116e-7,
    ];
    if x < 0.5 {
        std::f64::consts::PI.ln() - (std::f64::consts::PI * x).sin().ln() - lgamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let mut sum = 0.999_999_999_999_809_9;
        let t = x + 7.5;
        for (i, v) in p.iter().enumerate() {
            sum += v / (x + (i as f64) + 1.0);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + sum.ln()
    }
}

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
pub fn silhouette_score(partition: &Partition, dist_mat: &dyn DistanceMatrix) -> f64 {
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

/// Represents a set of coordinates for items: Item -> Vector
#[derive(Debug, Clone)]
pub struct Coordinates {
    pub data: HashMap<String, Vec<f64>>,
    pub dim: usize,
}

impl Coordinates {
    /// Load coordinates from a FeatureVector file.
    /// Format: Name <tab> Val1,Val2,Val3...
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
            let fv = FeatureVector::parse(&line);
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
pub fn davies_bouldin_score(partition: &Partition, coords: &Coordinates) -> f64 {
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
            let vec = coords.data.get(*item).unwrap();
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
            let vec = coords.data.get(*item).unwrap();
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
        let stat_i = stats.get(&i).unwrap();

        for &j in stats.keys() {
            if i == j {
                continue;
            }
            let stat_j = stats.get(&j).unwrap();

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
pub fn calinski_harabasz_score(partition: &Partition, coords: &Coordinates) -> f64 {
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
        let vec = coords.data.get(*item).unwrap();
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
            let vec = coords.data.get(*item).unwrap();
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
            let vec = coords.data.get(*item).unwrap();
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

/// Calculate Dunn Index
///
/// D = min(inter_cluster_dist) / max(intra_cluster_diameter)
///
/// Where:
/// - inter_cluster_dist: min distance between points in different clusters (Single Linkage)
/// - intra_cluster_diameter: max distance between points in same cluster (Complete Linkage / Diameter)
///
/// Higher values indicate better clustering.
pub fn dunn_score(partition: &Partition, dist_mat: &dyn DistanceMatrix) -> f64 {
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

        for j in i + 1..cluster_ids.len() {
            let cid_j = cluster_ids[j];
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
pub fn pbm_score(partition: &Partition, coords: &Coordinates) -> f64 {
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
        let vec = coords.data.get(*item).unwrap();
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
        let vec = coords.data.get(*item).unwrap();
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
            let vec = coords.data.get(*item).unwrap();
            for (d, val) in vec.iter().enumerate() {
                cluster_centroid[d] += val;
            }
        }
        for val in cluster_centroid.iter_mut() {
            *val /= n_members as f64;
        }
        
        // E_W contribution
        for item in members {
            let vec = coords.data.get(*item).unwrap();
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
pub fn ball_hall_score(partition: &Partition, coords: &Coordinates) -> f64 {
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
            let vec = coords.data.get(*item).unwrap();
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
            let vec = coords.data.get(*item).unwrap();
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
pub fn xie_beni_score(partition: &Partition, coords: &Coordinates) -> f64 {
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
            let vec = coords.data.get(*item).unwrap();
            for (d, val) in vec.iter().enumerate() {
                cluster_centroid[d] += val;
            }
        }
        for val in cluster_centroid.iter_mut() {
            *val /= n_members as f64;
        }
        
        // WGSS contribution
        for item in members {
            let vec = coords.data.get(*item).unwrap();
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
pub fn wemmert_gancarski_score(partition: &Partition, coords: &Coordinates) -> f64 {
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
            let vec = coords.data.get(*item).unwrap();
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
        
        let current_centroid = centroids.get(&cid).unwrap();
        let mut sum_r = 0.0;

        for item in members {
            let vec = coords.data.get(*item).unwrap();
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
pub fn c_index_score(partition: &Partition, dist_mat: &dyn DistanceMatrix) -> f64 {
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
pub fn gamma_score(partition: &Partition, dist_mat: &dyn DistanceMatrix) -> f64 {
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
            let is_diff = if partition[items[i]] != partition[items[j]] { 1.0 } else { 0.0 };
            
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
pub fn tau_score(partition: &Partition, dist_mat: &dyn DistanceMatrix) -> f64 {
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
        
        for k in i..j {
            if pairs[k].is_diff {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::phylo::tree::Tree;

    #[test]
    fn test_tree_distance() {
        // ((A:2,B:4)g:2,(C:2,((D:3,E:1)h:1,F:2)i:1)j:1)k;
        // Construct tree manually or parse from string?
        // We don't have direct access to parser here easily unless we export it.
        // Tree::from_newick is likely available if io module supports it.
        // Tree::to_newick is available.
        // Let's assume we can parse it.
        // But Tree::from_file is available.
        // Let's create a dummy tree manually.
        
        let mut tree = Tree::new();
        let root = tree.add_node(); // k
        tree.set_root(root);
        
        let g = tree.add_node();
        let j = tree.add_node();
        tree.add_child(root, g).unwrap();
        tree.add_child(root, j).unwrap();
        
        let a = tree.add_node();
        let b = tree.add_node();
        tree.get_node_mut(a).unwrap().name = Some("A".to_string());
        tree.get_node_mut(b).unwrap().name = Some("B".to_string());
        tree.add_child(g, a).unwrap();
        tree.add_child(g, b).unwrap();
        
        // Set lengths
        // A:2, B:4. g:2 (to k?)
        // If we treat edges as parent->child length.
        tree.get_node_mut(a).unwrap().length = Some(2.0);
        tree.get_node_mut(b).unwrap().length = Some(4.0);
        tree.get_node_mut(g).unwrap().length = Some(2.0);
        
        let td = TreeDistance::new(tree);
        
        // Distance A-B = 2+4 = 6.
        assert_eq!(td.get_distance("A", "B"), 6.0);
    }

    #[test]
    fn test_silhouette_score_simple() {
        // Data:
        // 0: 0.0 (C0)
        // 1: 1.0 (C1)
        // 2: 1.0 (C1)
        // 3: 2.0 (C1)
        // 4: 3.0 (C2)
        // 5: 3.0 (C2)

        let mut p = Partition::new();
        p.insert("0".to_string(), 0);
        p.insert("1".to_string(), 1);
        p.insert("2".to_string(), 1);
        p.insert("3".to_string(), 1);
        p.insert("4".to_string(), 2);
        p.insert("5".to_string(), 2);

        let names: Vec<String> = (0..6).map(|i| i.to_string()).collect();
        let mut dist_mat = NamedMatrix::new(names);
        let points: Vec<f32> = vec![0.0, 1.0, 1.0, 2.0, 3.0, 3.0];

        for i in 0..6 {
            for j in i + 1..6 {
                let d = (points[i] - points[j]).abs();
                let n1 = i.to_string();
                let n2 = j.to_string();
                dist_mat.set_by_name(&n1, &n2, d).unwrap();
            }
        }

        let score = silhouette_score(&p, &dist_mat);
        assert!((score - 0.5).abs() < 1e-6, "Score was {}", score);
    }

    #[test]
    fn test_silhouette_score_single_cluster() {
        let mut p = Partition::new();
        p.insert("0".to_string(), 0);
        p.insert("1".to_string(), 0);

        let names = vec!["0".to_string(), "1".to_string()];
        let mut dist_mat = NamedMatrix::new(names);
        dist_mat.set_by_name("0", "1", 1.0).unwrap();

        let score = silhouette_score(&p, &dist_mat);
        assert_eq!(score, 0.0);
    }

    #[test]
    fn test_silhouette_score_all_singletons() {
        // Sklearn behavior for all singletons is not strictly defined in docs but usually handled.
        // Our implementation returns 0.0 if n_clusters == n_samples
        let mut p = Partition::new();
        p.insert("0".to_string(), 0);
        p.insert("1".to_string(), 1);
        p.insert("2".to_string(), 2);

        let names = vec!["0".to_string(), "1".to_string(), "2".to_string()];
        let mut dist_mat = NamedMatrix::new(names);
        dist_mat.set_by_name("0", "1", 1.0).unwrap();
        dist_mat.set_by_name("0", "2", 1.0).unwrap();
        dist_mat.set_by_name("1", "2", 1.0).unwrap();

        let score = silhouette_score(&p, &dist_mat);
        assert_eq!(score, 0.0);
    }

    #[test]
    fn test_davies_bouldin_score_simple() {
        // Cluster 1: A(0,0), B(0,1) -> Centroid (0, 0.5), Scatter = 0.5
        // Cluster 2: C(5,0), D(5,1) -> Centroid (5, 0.5), Scatter = 0.5
        // M12 = 5.0
        // R12 = (0.5+0.5)/5.0 = 0.2
        // DB = (0.2 + 0.2)/2 = 0.2

        let mut p = Partition::new();
        p.insert("A".to_string(), 1);
        p.insert("B".to_string(), 1);
        p.insert("C".to_string(), 2);
        p.insert("D".to_string(), 2);

        let mut data = HashMap::new();
        data.insert("A".to_string(), vec![0.0, 0.0]);
        data.insert("B".to_string(), vec![0.0, 1.0]);
        data.insert("C".to_string(), vec![5.0, 0.0]);
        data.insert("D".to_string(), vec![5.0, 1.0]);

        let coords = Coordinates { data, dim: 2 };

        let score = davies_bouldin_score(&p, &coords);
        assert!((score - 0.2).abs() < 1e-6, "Score was {}", score);
    }

    #[test]
    fn test_evaluate_perfect() {
        let mut p1 = Partition::new();
        p1.insert("A".to_string(), 1);
        p1.insert("B".to_string(), 1);
        p1.insert("C".to_string(), 2);

        let mut p2 = Partition::new();
        p2.insert("A".to_string(), 10);
        p2.insert("B".to_string(), 10);
        p2.insert("C".to_string(), 20);

        let m = evaluate(&p1, &p2);
        assert_eq!(m.ari, 1.0);
        assert_eq!(m.ami, 1.0);
        assert_eq!(m.homogeneity, 1.0);
        assert_eq!(m.completeness, 1.0);
        assert_eq!(m.v_measure, 1.0);
        assert_eq!(m.fmi, 1.0);
        assert_eq!(m.nmi, 1.0);
    }

    #[test]
    fn test_evaluate_disjoint() {
        // P1: {A,B}, {C,D} -> Labels: 1, 1, 2, 2
        // P2: {A,C}, {B,D} -> Labels: 1, 2, 1, 2
        // Contingency table is uniform:
        //      P2_1(AC) P2_2(BD)
        // P1_1(AB)  1(A)     1(B)
        // P1_2(CD)  1(C)     1(D)
        //
        // This is perfectly independent (orthogonal).
        // MI = 0.0
        // NMI = 0.0
        // ARI = -0.5 (Worse than random?) Let's check calculation:
        // sum_nij_2 = 0
        // sum_a_2 = 1 + 1 = 2
        // sum_b_2 = 1 + 1 = 2
        // n_2 = binom(4, 2) = 6
        // E[Index] = (2 * 2) / 6 = 4/6 = 0.666
        // Max[Index] = (2 + 2) / 2 = 2
        // Index = 0
        // ARI = (0 - 0.666) / (2 - 0.666) = -0.666 / 1.333 = -0.5
        // FMI = TP / sqrt(2 * 2) = 0 / 2 = 0.0

        let mut p1 = Partition::new();
        p1.insert("A".to_string(), 1);
        p1.insert("B".to_string(), 1);
        p1.insert("C".to_string(), 2);
        p1.insert("D".to_string(), 2);

        let mut p2 = Partition::new();
        p2.insert("A".to_string(), 1);
        p2.insert("C".to_string(), 1);
        p2.insert("B".to_string(), 2);
        p2.insert("D".to_string(), 2);

        let m = evaluate(&p1, &p2);
        assert!((m.ari + 0.5).abs() < 1e-6);
        assert!(m.mi.abs() < 1e-6);
        assert!(m.nmi.abs() < 1e-6);
        assert!(m.fmi.abs() < 1e-6);
    }

    #[test]
    fn test_internal_indices_simple() {
        // Cluster 1: A(0,0), B(1,0) -> Centroid (0.5, 0)
        // Cluster 2: C(5,0), D(6,0) -> Centroid (5.5, 0)
        
        let mut p = Partition::new();
        p.insert("A".to_string(), 1);
        p.insert("B".to_string(), 1);
        p.insert("C".to_string(), 2);
        p.insert("D".to_string(), 2);

        let mut data = HashMap::new();
        data.insert("A".to_string(), vec![0.0, 0.0]);
        data.insert("B".to_string(), vec![1.0, 0.0]);
        data.insert("C".to_string(), vec![5.0, 0.0]);
        data.insert("D".to_string(), vec![6.0, 0.0]);

        let coords = Coordinates { data, dim: 2 };
        
        // Construct Distance Matrix for C-index
        let names = vec!["A".to_string(), "B".to_string(), "C".to_string(), "D".to_string()];
        let mut dist_mat = NamedMatrix::new(names);
        // Distances:
        // A-B: 1.0 (Intra)
        // C-D: 1.0 (Intra)
        // A-C: 5.0
        // A-D: 6.0
        // B-C: 4.0
        // B-D: 5.0
        // Sorted: 1, 1, 4, 5, 5, 6
        // N_W = 2 (A-B, C-D)
        // S_W = 1.0 + 1.0 = 2.0
        // S_min (sum of smallest 2) = 1.0 + 1.0 = 2.0
        // S_max (sum of largest 2) = 6.0 + 5.0 = 11.0
        // C-index = (2.0 - 2.0) / (11.0 - 2.0) = 0.0
        dist_mat.set_by_name("A", "B", 1.0).unwrap();
        dist_mat.set_by_name("C", "D", 1.0).unwrap();
        dist_mat.set_by_name("A", "C", 5.0).unwrap();
        dist_mat.set_by_name("A", "D", 6.0).unwrap();
        dist_mat.set_by_name("B", "C", 4.0).unwrap();
        dist_mat.set_by_name("B", "D", 5.0).unwrap();
        
        let c_index = c_index_score(&p, &dist_mat);
        assert_eq!(c_index, 0.0);

        // PBM:
        // Global Centroid: (12/4, 0) = (3, 0)
        // E_T: |0-3| + |1-3| + |5-3| + |6-3| = 3 + 2 + 2 + 3 = 10
        // E_W: 
        //   C1: |0-0.5| + |1-0.5| = 0.5 + 0.5 = 1.0
        //   C2: |5-5.5| + |6-5.5| = 0.5 + 0.5 = 1.0
        //   Total E_W = 2.0
        // D_B: |0.5 - 5.5| = 5.0
        // K = 2
        // PBM = ( 1/2 * 10 / 2.0 * 5.0 )^2 = ( 0.5 * 5 * 5 )^2 = (12.5)^2 = 156.25
        let pbm = pbm_score(&p, &coords);
        assert!((pbm - 156.25).abs() < 1e-6, "PBM was {}", pbm);
        
        // Ball-Hall:
        // C1 mean dispersion: (|0-0.5|^2 + |1-0.5|^2) / 2 = (0.25+0.25)/2 = 0.25
        // C2 mean dispersion: (|5-5.5|^2 + |6-5.5|^2) / 2 = (0.25+0.25)/2 = 0.25
        // BH = (0.25 + 0.25) / 2 = 0.25
        let bh = ball_hall_score(&p, &coords);
        assert!((bh - 0.25).abs() < 1e-6, "Ball-Hall was {}", bh);
        
        // Xie-Beni:
        // WGSS = (0.25+0.25) + (0.25+0.25) = 1.0
        // min_sq_dist = (5.0)^2 = 25.0
        // N = 4
        // XB = 1.0 / (4 * 25.0) = 0.01
        let xb = xie_beni_score(&p, &coords);
        assert!((xb - 0.01).abs() < 1e-6, "Xie-Beni was {}", xb);
        
        // Wemmert-Gancarski:
        // C1 centroid G1=(0.5,0), C2 centroid G2=(5.5,0)
        // A(0,0): ||A-G1||=0.5, ||A-G2||=5.5. R(A)=0.5/5.5 = 1/11
        // B(1,0): ||B-G1||=0.5, ||B-G2||=4.5. R(B)=0.5/4.5 = 1/9
        // C(5,0): ||C-G2||=0.5, ||C-G1||=4.5. R(C)=0.5/4.5 = 1/9
        // D(6,0): ||D-G2||=0.5, ||D-G1||=5.5. R(D)=0.5/5.5 = 1/11
        // Mean R(C1) = (1/11 + 1/9)/2 = (9/99 + 11/99)/2 = 20/99/2 = 10/99
        // Mean R(C2) = 10/99
        // J1 = 1 - 10/99 = 89/99
        // J2 = 1 - 10/99 = 89/99
        // J = (2/4)*J1 + (2/4)*J2 = 0.5*J1 + 0.5*J2 = 89/99 = 0.898989...
        let wg = wemmert_gancarski_score(&p, &coords);
        let expected_wg = 89.0 / 99.0;
        assert!((wg - expected_wg).abs() < 1e-6, "WG was {}, expected {}", wg, expected_wg);
    }
}
