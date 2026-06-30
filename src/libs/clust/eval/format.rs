//! Output formatting helpers for clustering evaluation metrics.

use super::coordinates::Coordinates;
use super::distance::DistanceMatrix;
use super::pairwise::Metrics;
use super::LabelMap;
use super::{
    ball_hall_score, c_index_score, calinski_harabasz_score, davies_bouldin_score, dunn_score,
    gamma_score, pbm_score, silhouette_score, tau_score, wemmert_gancarski_score, xie_beni_score,
};

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

/// Distance-based metric values, in DISTANCE_METRIC_NAMES order.
pub fn distance_metric_values(partition: &LabelMap, dist_mat: &dyn DistanceMatrix) -> Vec<f64> {
    vec![
        silhouette_score(partition, dist_mat),
        dunn_score(partition, dist_mat),
        c_index_score(partition, dist_mat),
        gamma_score(partition, dist_mat),
        tau_score(partition, dist_mat),
    ]
}

/// Coordinate-based metric values, in COORD_METRIC_NAMES order.
pub fn coord_metric_values(partition: &LabelMap, coords: &Coordinates) -> Vec<f64> {
    vec![
        davies_bouldin_score(partition, coords),
        calinski_harabasz_score(partition, coords),
        pbm_score(partition, coords),
        ball_hall_score(partition, coords),
        xie_beni_score(partition, coords),
        wemmert_gancarski_score(partition, coords),
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
