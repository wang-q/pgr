//! Topology report over the induced coarse GFA graph.

use super::dsu::Dsu;
use super::PafGraph;
use std::collections::HashMap;

/// Coarse GFA topology report.
#[derive(Debug)]
pub struct GraphReport {
    pub segments: usize,
    pub links: usize,
    pub paths: usize,
    pub path_steps: usize,
    pub total_segment_bp: usize,
    pub segment_len_min: usize,
    pub segment_len_mean: f64,
    pub segment_len_median: usize,
    pub segment_len_max: usize,
    pub node_coverage_mean: f64,
    pub node_coverage_median: usize,
    pub singleton_nodes: usize,
    pub reused_nodes: usize,
    pub reused_nodes_cross_path: usize,
    pub components: usize,
    pub largest_component_nodes: usize,
    pub tips: usize,
    pub isolated_nodes: usize,
    pub self_loop_edges: usize,
    pub path_len_steps_min: usize,
    pub path_len_steps_median: usize,
    pub path_len_steps_max: usize,
    pub path_len_bp_min: usize,
    pub path_len_bp_median: usize,
    pub path_len_bp_max: usize,
}

impl GraphReport {
    /// Write the report as TSV (key<TAB>value) to a writer.
    pub fn write_tsv<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        let fields: Vec<(&str, String)> = vec![
            ("segments", self.segments.to_string()),
            ("links", self.links.to_string()),
            ("paths", self.paths.to_string()),
            ("path_steps", self.path_steps.to_string()),
            ("total_segment_bp", self.total_segment_bp.to_string()),
            ("segment_len_min", self.segment_len_min.to_string()),
            ("segment_len_mean", format!("{:.2}", self.segment_len_mean)),
            ("segment_len_median", self.segment_len_median.to_string()),
            ("segment_len_max", self.segment_len_max.to_string()),
            (
                "node_coverage_mean",
                format!("{:.4}", self.node_coverage_mean),
            ),
            (
                "node_coverage_median",
                self.node_coverage_median.to_string(),
            ),
            ("singleton_nodes", self.singleton_nodes.to_string()),
            ("reused_nodes", self.reused_nodes.to_string()),
            (
                "reused_nodes_cross_path",
                self.reused_nodes_cross_path.to_string(),
            ),
            ("components", self.components.to_string()),
            (
                "largest_component_nodes",
                self.largest_component_nodes.to_string(),
            ),
            ("tips", self.tips.to_string()),
            ("isolated_nodes", self.isolated_nodes.to_string()),
            ("self_loop_edges", self.self_loop_edges.to_string()),
            ("path_len_steps_min", self.path_len_steps_min.to_string()),
            (
                "path_len_steps_median",
                self.path_len_steps_median.to_string(),
            ),
            ("path_len_steps_max", self.path_len_steps_max.to_string()),
            ("path_len_bp_min", self.path_len_bp_min.to_string()),
            ("path_len_bp_median", self.path_len_bp_median.to_string()),
            ("path_len_bp_max", self.path_len_bp_max.to_string()),
        ];
        for (k, v) in &fields {
            writeln!(w, "{k}\t{v}")?;
        }
        Ok(())
    }
}

impl PafGraph {
    /// Compute a topology report over the induced graph.
    pub fn report(&self) -> GraphReport {
        let segments = self.node_seqs.len();
        let links = self.edges.len();
        let paths = self.paths.len();
        let path_steps: usize = self.paths.iter().map(|(_, s)| s.len()).sum();
        let total_segment_bp: usize = self.node_lens.iter().sum();

        // Segment length distribution.
        let mut seg_lens: Vec<usize> = self.node_lens.clone();
        seg_lens.sort_unstable();
        let segment_len_min = seg_lens.first().copied().unwrap_or(0);
        let segment_len_max = seg_lens.last().copied().unwrap_or(0);
        let segment_len_mean = if segments > 0 {
            total_segment_bp as f64 / segments as f64
        } else {
            0.0
        };
        let segment_len_median = median_sorted(&seg_lens);

        // Per-node path-step coverage (how many path steps visit each node).
        let mut coverage = vec![0usize; segments];
        for (_, steps) in &self.paths {
            for st in steps {
                if (st.node as usize) < coverage.len() {
                    coverage[st.node as usize] += 1;
                }
            }
        }
        let node_coverage_mean = if segments > 0 {
            coverage.iter().sum::<usize>() as f64 / segments as f64
        } else {
            0.0
        };
        let singleton_nodes = coverage.iter().filter(|&&c| c == 1).count();
        let reused_nodes = coverage.iter().filter(|&&c| c > 1).count();

        // Per-path node set (to count cross-path reuse).
        let mut node_path_sets: Vec<std::collections::HashSet<u32>> = (0..segments)
            .map(|_| std::collections::HashSet::new())
            .collect();
        for (pi, (_, steps)) in self.paths.iter().enumerate() {
            for st in steps {
                if (st.node as usize) < node_path_sets.len() {
                    node_path_sets[st.node as usize].insert(pi as u32);
                }
            }
        }
        let reused_nodes_cross_path = node_path_sets.iter().filter(|s| s.len() > 1).count();

        let mut cov_sorted = coverage.clone();
        cov_sorted.sort_unstable();
        let node_coverage_median = median_sorted(&cov_sorted);

        // Node degree (undirected: count both endpoints; self-loop contributes 2).
        let mut degree = vec![0usize; segments];
        let mut self_loop_edges = 0usize;
        for e in &self.edges {
            let f = e.from as usize;
            let t = e.to as usize;
            if f == t {
                self_loop_edges += 1;
                degree[f] += 2;
            } else {
                if f < degree.len() {
                    degree[f] += 1;
                }
                if t < degree.len() {
                    degree[t] += 1;
                }
            }
        }
        let tips = degree.iter().filter(|&&d| d == 1).count();
        let isolated_nodes = degree.iter().filter(|&&d| d == 0).count();

        // Connected components (undirected DSU over edges).
        let mut dsu = Dsu::new(segments);
        for e in &self.edges {
            if (e.from as usize) < segments && (e.to as usize) < segments && e.from != e.to {
                dsu.union(e.from as usize, e.to as usize);
            }
        }
        let mut comp_size: HashMap<usize, usize> = HashMap::new();
        for i in 0..segments {
            *comp_size.entry(dsu.find(i)).or_insert(0) += 1;
        }
        let components = comp_size.len();
        let largest_component_nodes = comp_size.values().copied().max().unwrap_or(0);

        // Path length distribution (steps and bp).
        let mut path_steps_lens: Vec<usize> = self.paths.iter().map(|(_, s)| s.len()).collect();
        path_steps_lens.sort_unstable();
        let path_len_steps_min = path_steps_lens.first().copied().unwrap_or(0);
        let path_len_steps_max = path_steps_lens.last().copied().unwrap_or(0);
        let path_len_steps_median = median_sorted(&path_steps_lens);

        let mut path_bp_lens: Vec<usize> = self
            .paths
            .iter()
            .map(|(name, steps)| {
                let mut bp = 0usize;
                for st in steps {
                    if (st.node as usize) < self.node_seqs.len() {
                        bp += self.node_seqs[st.node as usize].len();
                    }
                }
                let _ = name;
                bp
            })
            .collect();
        path_bp_lens.sort_unstable();
        let path_len_bp_min = path_bp_lens.first().copied().unwrap_or(0);
        let path_len_bp_max = path_bp_lens.last().copied().unwrap_or(0);
        let path_len_bp_median = median_sorted(&path_bp_lens);

        GraphReport {
            segments,
            links,
            paths,
            path_steps,
            total_segment_bp,
            segment_len_min,
            segment_len_mean,
            segment_len_median,
            segment_len_max,
            node_coverage_mean,
            node_coverage_median,
            singleton_nodes,
            reused_nodes,
            reused_nodes_cross_path,
            components,
            largest_component_nodes,
            tips,
            isolated_nodes,
            self_loop_edges,
            path_len_steps_min,
            path_len_steps_median,
            path_len_steps_max,
            path_len_bp_min,
            path_len_bp_median,
            path_len_bp_max,
        }
    }
}

/// Median of a sorted slice (0 if empty).
fn median_sorted(sorted: &[usize]) -> usize {
    if sorted.is_empty() {
        return 0;
    }
    let n = sorted.len();
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2
    }
}
