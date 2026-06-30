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
pub use format::{
    coord_metric_values, distance_metric_values, external_metric_values, format_metrics_row,
    COORD_METRIC_NAMES, DISTANCE_METRIC_NAMES, EXTERNAL_METRIC_NAMES,
};
pub use pairwise::{evaluate, Metrics};
pub use partition::{load_batch_partitions, load_partition, remove_singletons, PartitionFormat};

use std::collections::HashMap;
use std::io::Write;

/// Represents a clustering partition: Item -> ClusterID
pub type LabelMap = HashMap<String, u32>;

/// Evaluation target (clap-free). The caller loads resources and constructs
/// this enum; `run_single`/`run_batch` consume it.
pub enum EvalTarget<'a> {
    /// External partition for ARI/AMI/V-Measure.
    External(&'a LabelMap),
    /// Distance matrix for silhouette/dunn/etc.
    Matrix(&'a dyn DistanceMatrix),
    /// Coordinates for davies-bouldin/etc.
    Coords(&'a Coordinates),
}

/// Run single evaluation of `p1` against `target`. Writes header + one row.
pub fn run_single<W: Write>(
    p1: &LabelMap,
    target: EvalTarget<'_>,
    writer: &mut W,
) -> anyhow::Result<()> {
    match target {
        EvalTarget::External(p2) => {
            let metrics = evaluate(p1, p2);
            writeln!(writer, "{}", EXTERNAL_METRIC_NAMES.join("\t"))?;
            writeln!(
                writer,
                "{}",
                format_metrics_row(&external_metric_values(&metrics))
            )?;
        }
        EvalTarget::Matrix(d) => {
            let values = distance_metric_values(p1, d);
            writeln!(writer, "{}", DISTANCE_METRIC_NAMES.join("\t"))?;
            writeln!(writer, "{}", format_metrics_row(&values))?;
        }
        EvalTarget::Coords(c) => {
            let values = coord_metric_values(p1, c);
            writeln!(writer, "{}", COORD_METRIC_NAMES.join("\t"))?;
            writeln!(writer, "{}", format_metrics_row(&values))?;
        }
    }
    Ok(())
}

/// Run batch evaluation. Writes a dynamic header (Group + each target's metric
/// names) followed by one row per batch group.
pub fn run_batch<W: Write>(
    batches: Vec<(String, LabelMap)>,
    targets: &[EvalTarget<'_>],
    writer: &mut W,
) -> anyhow::Result<()> {
    if targets.is_empty() {
        anyhow::bail!("at least one evaluation target required");
    }

    let mut header: Vec<&str> = vec!["Group"];
    for t in targets {
        match t {
            EvalTarget::External(_) => header.extend_from_slice(EXTERNAL_METRIC_NAMES),
            EvalTarget::Matrix(_) => header.extend_from_slice(DISTANCE_METRIC_NAMES),
            EvalTarget::Coords(_) => header.extend_from_slice(COORD_METRIC_NAMES),
        }
    }
    writeln!(writer, "{}", header.join("\t"))?;

    for (group, p1) in batches {
        let mut row: Vec<String> = vec![group];
        for t in targets {
            match t {
                EvalTarget::External(p2) => {
                    let metrics = evaluate(&p1, p2);
                    row.push(format_metrics_row(&external_metric_values(&metrics)));
                }
                EvalTarget::Matrix(d) => {
                    let values = distance_metric_values(&p1, *d);
                    row.push(format_metrics_row(&values));
                }
                EvalTarget::Coords(c) => {
                    let values = coord_metric_values(&p1, c);
                    row.push(format_metrics_row(&values));
                }
            }
        }
        writeln!(writer, "{}", row.join("\t"))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests;
