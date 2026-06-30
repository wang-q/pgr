//! Coarse GFA graph induction from PAF (seqwish-style segment-level DSU).
//!
//! Splits each alignment at indels >= `min_var_len` into "match segments",
//! unions aligned segments via a disjoint-set union (transitive closure),
//! derives graph nodes (DSU classes) + edges (path adjacencies) + novel
//! segments (unaligned gaps), and emits GFA v1.0 (S/L/P).

mod dsu;
mod gfa;
mod segment;
#[cfg(test)]
mod tests;

use super::cigar::extract_cigar;
use super::parser::parse_paf;
use dsu::Dsu;
use segment::{id_to_name_local, novel_node_for, split_alignment, AlignmentLink, Segment};
use std::collections::HashMap;
use std::io::BufRead;

/// A graph edge (GFA L line).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Edge {
    pub from: u32,
    pub from_orient: char,
    pub to: u32,
    pub to_orient: char,
}

/// A path entry (one node visit in a GFA P line).
#[derive(Clone, Debug)]
pub struct PathStep {
    pub node: u32,
    pub orient: char,
}

/// Induced coarse graph, ready for GFA emission.
pub struct PafGraph {
    /// Node sequences (one per DSU class), indexed by node id.
    /// Empty when built without FASTA (topology-only mode).
    pub node_seqs: Vec<Vec<u8>>,
    /// Node lengths (bp). Populated even in topology-only mode from segment coords.
    pub node_lens: Vec<usize>,
    /// Per-node rGFA origin: (source sequence name, 0-based start offset).
    pub node_origins: Vec<(String, i32)>,
    /// Edges (deduplicated).
    pub edges: Vec<Edge>,
    /// Paths: (sequence name, steps).
    pub paths: Vec<(String, Vec<PathStep>)>,
}

impl PafGraph {
    /// Build a coarse GFA graph from a PAF reader + per-sequence FASTA bytes.
    ///
    /// `seqs` maps sequence name -> (forward-strand bytes). `min_var_len` is the
    /// minimum indel length to split at (smaller indels stay within a segment).
    pub fn build<R: BufRead>(
        paf_reader: R,
        seqs: Option<&HashMap<String, Vec<u8>>>,
        min_var_len: i32,
    ) -> anyhow::Result<Self> {
        let records = parse_paf(paf_reader)?;

        // Build name -> id map (reuse order from records; fall back to seqs keys).
        let mut name_to_id: HashMap<String, u32> = HashMap::new();
        let register = |name: &str, map: &mut HashMap<String, u32>| -> u32 {
            if let Some(&id) = map.get(name) {
                id
            } else {
                let id = map.len() as u32;
                map.insert(name.to_string(), id);
                id
            }
        };
        for r in &records {
            register(&r.target_name, &mut name_to_id);
            register(&r.query_name, &mut name_to_id);
        }
        if let Some(seqs_map) = seqs {
            for name in seqs_map.keys() {
                register(name, &mut name_to_id);
            }
        }

        // ── Stage 1: split alignments at SV breakpoints → segments + links ──
        let mut segments: Vec<Segment> = Vec::new();
        let mut links: Vec<AlignmentLink> = Vec::new();

        for rec in &records {
            let tid = name_to_id[&rec.target_name];
            let qid = name_to_id[&rec.query_name];
            let reverse = rec.strand == '-';
            let cigar = extract_cigar(&rec.tags)?;
            split_alignment(
                tid,
                rec.target_start as i32,
                rec.target_end as i32,
                qid,
                rec.query_start as i32,
                rec.query_end as i32,
                rec.query_length as i32,
                reverse,
                &cigar,
                min_var_len,
                &mut segments,
                &mut links,
            );
        }

        // ── Stage 2: DSU transitive closure ──
        let mut dsu = Dsu::new(segments.len());
        for link in &links {
            dsu.union(link.a, link.b);
        }

        // Assign node ids by DSU root (sorted by min segment start for stability).
        let mut root_to_node: HashMap<usize, u32> = HashMap::new();
        // Collect (root, min_seq_id, min_start, seg_idx) to order nodes.
        let mut root_info: Vec<(usize, u32, i32, usize)> = Vec::new();
        for (i, seg) in segments.iter().enumerate() {
            let root = dsu.find(i);
            root_info.push((root, seg.seq_id, seg.start, i));
        }
        // Sort by (seq_id, start) so node 0 is the earliest segment — stable GFA.
        root_info.sort_by_key(|&(root, sid, start, _)| (root, sid, start));
        for &(root, _, _, _) in &root_info {
            let node_id = root_to_node.len() as u32;
            root_to_node.entry(root).or_insert(node_id);
        }

        let num_nodes = root_to_node.len() as u32;
        // Map segment index -> node id.
        let seg_node: Vec<u32> = (0..segments.len())
            .map(|i| root_to_node[&dsu.find(i)])
            .collect();

        // ── Stage 3: node sequences (first-seen segment's forward strand) ──
        let mut node_seqs: Vec<Vec<u8>> = vec![Vec::new(); num_nodes as usize];
        let mut node_origins: Vec<(String, i32)> = vec![(String::new(), 0); num_nodes as usize];
        let mut node_filled: Vec<bool> = vec![false; num_nodes as usize];
        // Walk segments in sorted order (same sort as node assignment) for stability.
        for &(_, _, _, seg_idx) in &root_info {
            let node = seg_node[seg_idx] as usize;
            if node_filled[node] {
                continue;
            }
            let seg = &segments[seg_idx];
            if let Some(name) = id_to_name_local(&name_to_id, seg.seq_id) {
                if let Some(seq_bytes) = seqs.and_then(|m| m.get(name)) {
                    let s = seg.start.max(0) as usize;
                    let e = (seg.end as usize).min(seq_bytes.len());
                    if s < e {
                        node_seqs[node] = seq_bytes[s..e].to_vec();
                        node_origins[node] = (name.to_string(), seg.start);
                        node_filled[node] = true;
                    }
                }
            }
        }

        // ── Stage 4: path construction per sequence ──
        let num_seqs = name_to_id.len() as u32;
        let mut paths: Vec<(String, Vec<PathStep>)> = Vec::new();
        let mut edges: Vec<Edge> = Vec::new();

        for sid in 0..num_seqs {
            let name = id_to_name_local(&name_to_id, sid)
                .unwrap_or("?")
                .to_string();
            // Collect segments on this sequence, sorted by start.
            let mut segs_on_seq: Vec<(usize, &Segment)> = segments
                .iter()
                .enumerate()
                .filter(|(_, s)| s.seq_id == sid)
                .collect();
            segs_on_seq.sort_by_key(|(_, s)| s.start);

            let seq_len = seqs
                .and_then(|m| m.get(&name))
                .map(|v| v.len())
                .unwrap_or(0) as i32;
            let has_seqs = seqs.is_some();

            let mut steps: Vec<PathStep> = Vec::new();
            let mut cursor = 0i32;

            for &(seg_idx, seg) in &segs_on_seq {
                // Novel segment for the gap before this aligned segment.
                if has_seqs && seg.start > cursor {
                    let novel_node = novel_node_for(
                        &mut node_seqs,
                        &mut node_origins,
                        sid,
                        cursor,
                        seg.start,
                        seqs,
                        &name_to_id,
                    );
                    steps.push(PathStep {
                        node: novel_node,
                        orient: '+',
                    });
                }
                // The aligned segment's node.
                let node = seg_node[seg_idx];
                // Orientation: segments are always forward coords; the alignment
                // orientation affects the *edge*, not the path orientation here.
                // For simplicity, path orientation is '+' (forward traversal).
                // (Reverse-strand alignments produce forward-coord segments;
                //  a full rGFA would tag orientation, but coarse GFA uses '+'.)
                steps.push(PathStep { node, orient: '+' });
                cursor = seg.end.max(cursor);
            }
            // Trailing novel segment.
            if has_seqs && cursor < seq_len {
                let novel_node = novel_node_for(
                    &mut node_seqs,
                    &mut node_origins,
                    sid,
                    cursor,
                    seq_len,
                    seqs,
                    &name_to_id,
                );
                steps.push(PathStep {
                    node: novel_node,
                    orient: '+',
                });
            }

            // Derive edges from consecutive steps.
            for w in steps.windows(2) {
                let e = Edge {
                    from: w[0].node,
                    from_orient: w[0].orient,
                    to: w[1].node,
                    to_orient: w[1].orient,
                };
                if !edges.contains(&e) {
                    edges.push(e);
                }
            }

            if !steps.is_empty() {
                paths.push((name, steps));
            } else if has_seqs && seq_len > 0 {
                // No alignments at all: whole sequence is one novel node.
                let novel_node = novel_node_for(
                    &mut node_seqs,
                    &mut node_origins,
                    sid,
                    0,
                    seq_len,
                    seqs,
                    &name_to_id,
                );
                paths.push((
                    name,
                    vec![PathStep {
                        node: novel_node,
                        orient: '+',
                    }],
                ));
            }
        }

        // Compute node lengths from segment coords (works even without FASTA).
        let mut node_lens: Vec<usize> = vec![0; num_nodes as usize];
        for &(_, _, _, seg_idx) in &root_info {
            let node = seg_node[seg_idx] as usize;
            if node_lens[node] == 0 {
                let seg = &segments[seg_idx];
                node_lens[node] = (seg.end - seg.start).max(0) as usize;
            }
        }
        // Extend for novel nodes added during path construction.
        while node_lens.len() < node_seqs.len() {
            let i = node_lens.len();
            let len = node_seqs.get(i).map(|s| s.len()).unwrap_or(0);
            node_lens.push(len);
        }
        // Override with actual sequence lengths when available.
        for (i, seq) in node_seqs.iter().enumerate() {
            if !seq.is_empty() {
                node_lens[i] = seq.len();
            }
        }

        Ok(PafGraph {
            node_seqs,
            node_lens,
            node_origins,
            edges,
            paths,
        })
    }

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
