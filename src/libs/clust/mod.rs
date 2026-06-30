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
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 {
            continue;
        }
        let a = names.insert_full(fields[0].to_string()).0;
        let b = names.insert_full(fields[1].to_string()).0;
        graph.add_edge(a, b, ());
    }

    let scc = petgraph::algo::tarjan_scc(&graph);
    let names_vec: Vec<String> = names.iter().cloned().collect();
    Ok((names_vec, scc))
}
