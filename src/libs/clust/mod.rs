pub mod dbscan;
pub mod eval;
pub mod feature;
pub mod format;
pub mod hier;
pub mod k_medoids;
pub mod mcl;
pub mod medoid;
pub mod nj;
pub mod tree_cut;
pub mod upgma;

use anyhow::Result;
use indexmap::IndexSet;
use std::io::BufRead;

/// Load pairwise relations from a TSV reader and compute connected components.
///
/// Returns `(names, components)` where `names[i]` is the i-th node's name and
/// `components` is a Vec of Vecs of node indices (one Vec per component).
pub fn connected_components<R: BufRead>(reader: R) -> Result<(Vec<String>, Vec<Vec<usize>>)> {
    let mut names = IndexSet::new();
    let mut graph = petgraph::graphmap::UnGraphMap::<_, ()>::new();

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 {
            anyhow::bail!(
                "invalid pairwise relation line (expected at least 2 tab-separated columns): {}",
                line
            );
        }
        let a = names.insert_full(fields[0].to_string()).0;
        let b = names.insert_full(fields[1].to_string()).0;
        graph.add_edge(a, b, ());
    }

    let scc = petgraph::algo::tarjan_scc(&graph);
    let names_vec: Vec<String> = names.iter().cloned().collect();
    Ok((names_vec, scc))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    fn sort_components(components: &mut Vec<Vec<usize>>) {
        for c in components.iter_mut() {
            c.sort();
        }
        components.sort();
    }

    #[test]
    fn test_connected_components_empty() -> anyhow::Result<()> {
        let data = "";
        let reader = BufReader::new(data.as_bytes());
        let (names, components) = connected_components(reader)?;
        assert!(names.is_empty());
        assert!(components.is_empty());
        Ok(())
    }

    #[test]
    fn test_connected_components_single_pair() -> anyhow::Result<()> {
        let data = "A\tB\n";
        let reader = BufReader::new(data.as_bytes());
        let (names, mut components) = connected_components(reader)?;
        assert_eq!(names, vec!["A", "B"]);
        sort_components(&mut components);
        assert_eq!(components, vec![vec![0, 1]]);
        Ok(())
    }

    #[test]
    fn test_connected_components_multiple() -> anyhow::Result<()> {
        // Two disconnected components: {A,B} and {C,D,E}
        let data = "A\tB\nC\tD\nD\tE\n";
        let reader = BufReader::new(data.as_bytes());
        let (names, mut components) = connected_components(reader)?;
        assert_eq!(names, vec!["A", "B", "C", "D", "E"]);
        sort_components(&mut components);
        assert_eq!(components, vec![vec![0, 1], vec![2, 3, 4]]);
        Ok(())
    }

    #[test]
    fn test_connected_components_malformed() {
        let data = "A\n";
        let reader = BufReader::new(data.as_bytes());
        let result = connected_components(reader);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(
            err.contains("pairwise relation"),
            "unexpected error: {}",
            err
        );
    }
}
