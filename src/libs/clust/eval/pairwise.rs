use super::LabelMap;
use std::collections::{HashMap, HashSet};

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

pub fn evaluate(p1: &LabelMap, p2: &LabelMap) -> Metrics {
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

fn normalize_labels(p: &LabelMap, keys: &[&String]) -> (Vec<u32>, usize) {
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
        -176.615_029_162_140_6,
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
