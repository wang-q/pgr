//! Taxonomy TSV parsing for tree condensation pipelines.

use std::collections::{BTreeMap, BTreeSet};
use std::io::BufRead;

use super::parser::newick_safe;

/// Taxonomy table: per-node terms grouped by rank, plus unique groups per rank.
#[derive(Debug, Default)]
pub struct TaxonomyTable {
    /// node_name -> terms per rank (None if column missing or empty).
    pub taxon_map: BTreeMap<String, Vec<Option<String>>>,
    /// Unique, sorted, NA-filtered group names per rank.
    pub all_groups: Vec<Vec<String>>,
}

/// Read taxonomy TSV from `reader` and filter to `leaf_names`.
///
/// `ranks` is a list of 1-based column indices into the TSV. Lines with fewer
/// than 2 columns are skipped with a warning. Nodes not in `leaf_names` are
/// ignored. Group lists are deduplicated, sorted, and filtered to drop `"NA"`.
pub fn read_taxonomy<R: BufRead>(
    reader: R,
    ranks: &[usize],
    leaf_names: &BTreeSet<String>,
) -> anyhow::Result<TaxonomyTable> {
    let mut taxon_map: BTreeMap<String, Vec<Option<String>>> = BTreeMap::new();
    let mut all_groups: Vec<Vec<String>> = vec![vec![]; ranks.len()];

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            log::warn!("skipping line with <2 columns: {}", line);
            continue;
        }
        let node_name = parts[0].to_string();
        if !leaf_names.contains(&node_name) {
            continue;
        }
        let mut terms: Vec<Option<String>> = Vec::with_capacity(ranks.len());

        for (i, rank_col) in ranks.iter().enumerate() {
            let rank_idx = rank_col.saturating_sub(1);
            let term = parts.get(rank_idx).map(|s| newick_safe(s));
            if let Some(t) = &term {
                all_groups[i].push(t.clone());
            }
            terms.push(term);
        }

        if terms.iter().any(|t| t.is_some()) {
            taxon_map.insert(node_name, terms);
        }
    }

    // Deduplicate, sort, and filter NA per rank.
    for groups in &mut all_groups {
        groups.sort();
        groups.dedup();
        groups.retain(|s| s != "NA");
    }

    Ok(TaxonomyTable {
        taxon_map,
        all_groups,
    })
}
