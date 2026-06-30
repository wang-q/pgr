//! Output formatting helpers for clustering evaluation metrics.

use super::pairwise::Metrics;

/// External (pairwise) evaluation metric names, in output column order.
pub const EXTERNAL_METRIC_NAMES: &[&str] = &[
    "ari",
    "ami",
    "homogeneity",
    "completeness",
    "v_measure",
    "fmi",
    "nmi",
    "mi",
    "ri",
    "jaccard",
    "precision",
    "recall",
];

/// Distance-based evaluation metric names, in output column order.
pub const DISTANCE_METRIC_NAMES: &[&str] = &["silhouette", "dunn", "c_index", "gamma", "tau"];

/// Coordinate-based evaluation metric names, in output column order.
pub const COORD_METRIC_NAMES: &[&str] = &[
    "davies_bouldin",
    "calinski_harabasz",
    "pbm",
    "ball_hall",
    "xie_beni",
    "wemmert_gancarski",
];

/// External metric values from a Metrics struct, in EXTERNAL_METRIC_NAMES order.
pub fn external_metric_values(m: &Metrics) -> Vec<f64> {
    vec![
        m.ari,
        m.ami,
        m.homogeneity,
        m.completeness,
        m.v_measure,
        m.fmi,
        m.nmi,
        m.mi,
        m.ri,
        m.jaccard,
        m.precision,
        m.recall,
    ]
}

/// Format a slice of f64 values as tab-separated `{:.6}` strings.
pub fn format_metrics_row(values: &[f64]) -> String {
    values
        .iter()
        .map(|v| format!("{:.6}", v))
        .collect::<Vec<_>>()
        .join("\t")
}
