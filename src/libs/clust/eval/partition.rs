use super::LabelMap;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Supported partition file formats for clustering evaluation.
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
///    - If 2 columns: ClusterID `tab` Item
///    - If > 2 columns: Treated as Cluster-based.
/// 3. Long-based: Treated as Batch LabelMap (returns empty map here, use load_batch_partitions).
pub fn load_partition<P: AsRef<Path>>(
    path: P,
    format: PartitionFormat,
) -> anyhow::Result<LabelMap> {
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

fn parse_pair_format(lines: &[String]) -> anyhow::Result<LabelMap> {
    let mut partition = HashMap::new();
    // We need to map string labels to u32 IDs
    let mut label_map: HashMap<String, u32> = HashMap::new();
    let mut next_id = 0;

    for line in lines {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 {
            anyhow::bail!("invalid pair format line (expected 2 columns): {}", line);
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

fn parse_cluster_format(lines: &[String]) -> anyhow::Result<LabelMap> {
    let mut partition = HashMap::new();
    let mut cluster_id = 0;

    for line in lines {
        let parts: Vec<&str> = line.split_whitespace().collect();
        cluster_id += 1;
        for item in parts {
            partition.insert(item.to_string(), cluster_id);
        }
    }
    Ok(partition)
}

/// Load batch partitions from a file in Long format.
/// Format: GroupID `tab` ClusterID `tab` SampleID
/// Returns a list of (GroupID, LabelMap).
pub fn load_batch_partitions<P: AsRef<Path>>(path: P) -> anyhow::Result<Vec<(String, LabelMap)>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut groups: Vec<String> = Vec::new();
    let mut group_indices: HashMap<String, usize> = HashMap::new();
    let mut partitions: Vec<LabelMap> = Vec::new();

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
            anyhow::bail!(
                "invalid long format line (expected 3 tab-separated columns): {}",
                line
            );
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

    let result = groups.into_iter().zip(partitions).collect();
    Ok(result)
}

/// Remove singleton clusters (clusters with only one member) from a partition.
pub fn remove_singletons(partition: &mut LabelMap) {
    let mut counts = HashMap::new();
    for cid in partition.values() {
        *counts.entry(*cid).or_insert(0) += 1;
    }
    partition.retain(|_, cid| counts[cid] > 1);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_partition_format_parse() {
        assert_eq!(
            "cluster".parse::<PartitionFormat>().unwrap(),
            PartitionFormat::Cluster
        );
        assert_eq!(
            "pair".parse::<PartitionFormat>().unwrap(),
            PartitionFormat::Pair
        );
        assert_eq!(
            "long".parse::<PartitionFormat>().unwrap(),
            PartitionFormat::Long
        );
        assert!("unknown".parse::<PartitionFormat>().is_err());
    }

    #[test]
    fn test_load_partition_cluster() -> anyhow::Result<()> {
        let mut file = tempfile::NamedTempFile::new()?;
        writeln!(file, "A B")?;
        writeln!(file, "C")?;
        writeln!(file, "# comment")?;

        let partition = load_partition(file.path(), PartitionFormat::Cluster)?;
        assert_eq!(partition.get("A"), Some(&1));
        assert_eq!(partition.get("B"), Some(&1));
        assert_eq!(partition.get("C"), Some(&2));
        assert_eq!(partition.len(), 3);
        Ok(())
    }

    #[test]
    fn test_load_partition_pair() -> anyhow::Result<()> {
        let mut file = tempfile::NamedTempFile::new()?;
        writeln!(file, "1\tA")?;
        writeln!(file, "1\tB")?;
        writeln!(file, "2\tC")?;

        let partition = load_partition(file.path(), PartitionFormat::Pair)?;
        assert_eq!(partition.get("A"), Some(&1));
        assert_eq!(partition.get("B"), Some(&1));
        assert_eq!(partition.get("C"), Some(&2));
        assert_eq!(partition.len(), 3);
        Ok(())
    }

    #[test]
    fn test_load_partition_pair_malformed() -> anyhow::Result<()> {
        let mut file = tempfile::NamedTempFile::new()?;
        writeln!(file, "A")?;
        writeln!(file, "B")?;

        let result = load_partition(file.path(), PartitionFormat::Pair);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("pair format"), "unexpected error: {}", err);
        Ok(())
    }

    #[test]
    fn test_load_partition_long_rejected() -> anyhow::Result<()> {
        let mut file = tempfile::NamedTempFile::new()?;
        writeln!(file, "g\t1\tA")?;

        let result = load_partition(file.path(), PartitionFormat::Long);
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn test_load_batch_partitions() -> anyhow::Result<()> {
        let mut file = tempfile::NamedTempFile::new()?;
        writeln!(file, "Group\tClusterID\tSampleID")?;
        writeln!(file, "g1\t1\tA")?;
        writeln!(file, "g1\t1\tB")?;
        writeln!(file, "g2\t1\tC")?;

        let batches = load_batch_partitions(file.path())?;
        assert_eq!(batches.len(), 2);

        let (g1, p1) = &batches[0];
        assert_eq!(g1, "g1");
        assert_eq!(p1.get("A"), Some(&1));
        assert_eq!(p1.get("B"), Some(&1));

        let (g2, p2) = &batches[1];
        assert_eq!(g2, "g2");
        assert_eq!(p2.get("C"), Some(&1));
        Ok(())
    }

    #[test]
    fn test_load_batch_partitions_malformed() -> anyhow::Result<()> {
        let mut file = tempfile::NamedTempFile::new()?;
        writeln!(file, "g1\t1\tA")?;
        writeln!(file, "g1\tA")?;

        let result = load_batch_partitions(file.path());
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("long format"), "unexpected error: {}", err);
        Ok(())
    }

    #[test]
    fn test_remove_singletons() {
        let mut partition = HashMap::new();
        partition.insert("A".to_string(), 1);
        partition.insert("B".to_string(), 1);
        partition.insert("C".to_string(), 2);

        remove_singletons(&mut partition);
        assert!(!partition.contains_key("C"));
        assert_eq!(partition.get("A"), Some(&1));
        assert_eq!(partition.get("B"), Some(&1));
    }
}
