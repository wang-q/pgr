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
}

impl std::str::FromStr for PartitionFormat {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "cluster" => Ok(PartitionFormat::Cluster),
            "pair" => Ok(PartitionFormat::Pair),
            _ => Err(format!("Unknown format: {}", s)),
        }
    }
}

/// Load a partition from a file.
/// Supports two formats:
/// 1. Cluster-based: Each line is a cluster, items separated by whitespace.
///    The first item is treated as the cluster representative/ID.
/// 2. Pair-based: Two columns.
///    - If 2 columns: Item <tab> ClusterID
///    - If > 2 columns: Treated as Cluster-based.
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

#[derive(Debug, Default)]
pub struct Metrics {
    pub ari: f64,
    pub ami: f64,
    pub homogeneity: f64,
    pub completeness: f64,
    pub v_measure: f64,
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

    // 4. Calculate Metrics
    let ari = calculate_ari(&table, &a_counts, &b_counts, n);
    let (homogeneity, completeness, v_measure) =
        calculate_v_measure(&table, &a_counts, &b_counts, n);
    let ami = calculate_ami(&table, &a_counts, &b_counts, n);

    Metrics {
        ari,
        ami,
        homogeneity,
        completeness,
        v_measure,
    }
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

fn calculate_ari(
    table: &HashMap<(u32, u32), usize>,
    a_counts: &[usize],
    b_counts: &[usize],
    n: usize,
) -> f64 {
    if a_counts.len() <= 1 || b_counts.len() <= 1 {
        // Special case: If both partitions have only 1 cluster, it's a perfect match (ARI=1.0).
        // If only one partition has 1 cluster, ARI is 0.0.
        if a_counts.len() == 1 && b_counts.len() == 1 {
            return 1.0;
        }
        return 0.0;
    }

    fn binom2(x: usize) -> f64 {
        if x < 2 {
            0.0
        } else {
            (x as f64 * (x as f64 - 1.0)) / 2.0
        }
    }

    let sum_nij_2: f64 = table.values().map(|&count| binom2(count)).sum();
    let sum_a_2: f64 = a_counts.iter().map(|&count| binom2(count)).sum();
    let sum_b_2: f64 = b_counts.iter().map(|&count| binom2(count)).sum();
    let n_2 = binom2(n);

    let expected_index = (sum_a_2 * sum_b_2) / n_2;
    let max_index = (sum_a_2 + sum_b_2) / 2.0;
    let index = sum_nij_2;

    if max_index - expected_index == 0.0 {
        return 0.0; // Avoid division by zero
    }

    (index - expected_index) / (max_index - expected_index)
}

fn calculate_v_measure(
    table: &HashMap<(u32, u32), usize>,
    a_counts: &[usize],
    b_counts: &[usize],
    n: usize,
) -> (f64, f64, f64) {
    if n == 0 {
        return (0.0, 0.0, 0.0);
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

    // Homogeneity = 1 - H(U|V) / H(U) = MI(U,V) / H(U)
    // Completeness = 1 - H(V|U) / H(V) = MI(U,V) / H(V)

    let homogeneity = if h_a == 0.0 { 1.0 } else { mi / h_a };
    let completeness = if h_b == 0.0 { 1.0 } else { mi / h_b };

    let v_measure = if homogeneity + completeness == 0.0 {
        0.0
    } else {
        2.0 * homogeneity * completeness / (homogeneity + completeness)
    };

    (homogeneity, completeness, v_measure)
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
