//! Coarse GFA graph induction from PAF (seqwish-style segment-level DSU).
//!
//! Splits each alignment at indels >= `min_var_len` into "match segments",
//! unions aligned segments via a disjoint-set union (transitive closure),
//! derives graph nodes (DSU classes) + edges (path adjacencies) + novel
//! segments (unaligned gaps), and emits GFA v1.0 (S/L/P).

use super::cigar::CigarOp;
use super::parser::parse_paf;
use std::collections::HashMap;
use std::io::BufRead;

/// A forward-orientation region on a sequence.
#[derive(Clone, Debug)]
struct Segment {
    seq_id: u32,
    start: i32, // 0-based, inclusive
    end: i32,   // exclusive
}

/// A bidirectional alignment link between two segments.
struct AlignmentLink {
    a: usize, // index into `segments`
    b: usize, // index into `segments`
    #[allow(dead_code)]
    reverse: bool, // true if b is reverse-complement of a (for future rGFA orientation)
}

/// Simple union-find (union-by-rank + path compression).
struct Dsu {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl Dsu {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }
    fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }
    fn union(&mut self, a: usize, b: usize) {
        let (mut ra, mut rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        if self.rank[ra] < self.rank[rb] {
            std::mem::swap(&mut ra, &mut rb);
        }
        self.parent[rb] = ra;
        if self.rank[ra] == self.rank[rb] {
            self.rank[ra] += 1;
        }
    }
}

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
    pub node_seqs: Vec<Vec<u8>>,
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
        seqs: &HashMap<String, Vec<u8>>,
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
        for name in seqs.keys() {
            register(name, &mut name_to_id);
        }

        // ── Stage 1: split alignments at SV breakpoints → segments + links ──
        let mut segments: Vec<Segment> = Vec::new();
        let mut links: Vec<AlignmentLink> = Vec::new();

        for rec in &records {
            let tid = name_to_id[&rec.target_name];
            let qid = name_to_id[&rec.query_name];
            let reverse = rec.strand == '-';
            let cigar = extract_cigar(&rec.tags);
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
        let mut node_filled: Vec<bool> = vec![false; num_nodes as usize];
        // Walk segments in sorted order (same sort as node assignment) for stability.
        for &(_, _, _, seg_idx) in &root_info {
            let node = seg_node[seg_idx] as usize;
            if node_filled[node] {
                continue;
            }
            let seg = &segments[seg_idx];
            if let Some(name) = id_to_name_local(&name_to_id, seg.seq_id) {
                if let Some(seq_bytes) = seqs.get(name) {
                    let s = seg.start.max(0) as usize;
                    let e = (seg.end as usize).min(seq_bytes.len());
                    if s < e {
                        node_seqs[node] = seq_bytes[s..e].to_vec();
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

            let seq_len = seqs.get(&name).map(|v| v.len()).unwrap_or(0) as i32;

            let mut steps: Vec<PathStep> = Vec::new();
            let mut cursor = 0i32;

            for &(seg_idx, seg) in &segs_on_seq {
                // Novel segment for the gap before this aligned segment.
                if seg.start > cursor {
                    let novel_node =
                        novel_node_for(&mut node_seqs, sid, cursor, seg.start, seqs, &name_to_id);
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
            if cursor < seq_len {
                let novel_node =
                    novel_node_for(&mut node_seqs, sid, cursor, seq_len, seqs, &name_to_id);
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
            } else if seq_len > 0 {
                // No alignments at all: whole sequence is one novel node.
                let novel_node = novel_node_for(&mut node_seqs, sid, 0, seq_len, seqs, &name_to_id);
                paths.push((
                    name,
                    vec![PathStep {
                        node: novel_node,
                        orient: '+',
                    }],
                ));
            }
        }

        Ok(PafGraph {
            node_seqs,
            edges,
            paths,
        })
    }

    /// Write GFA v1.0 to a writer (S + L + P lines).
    pub fn write_gfa<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        // S lines (1-based node ids in GFA convention).
        for (i, seq) in self.node_seqs.iter().enumerate() {
            let id = (i + 1) as u32;
            let s = String::from_utf8_lossy(seq);
            writeln!(w, "S\t{id}\t{s}")?;
        }
        // L lines.
        for e in &self.edges {
            writeln!(
                w,
                "L\t{}\t{}\t{}\t{}\t0M",
                e.from + 1,
                e.from_orient,
                e.to + 1,
                e.to_orient
            )?;
        }
        // P lines.
        for (name, steps) in &self.paths {
            let path_str: Vec<String> = steps
                .iter()
                .map(|s| format!("{}{}", s.node + 1, s.orient))
                .collect();
            let overlaps = vec!["0M"; steps.len().saturating_sub(1)];
            writeln!(
                w,
                "P\t{name}\t{}\t{}",
                path_str.join(","),
                overlaps.join(",")
            )?;
        }
        Ok(())
    }

    /// Compute a topology report over the induced graph (V6 graph quality).
    pub fn report(&self) -> GraphReport {
        let segments = self.node_seqs.len();
        let links = self.edges.len();
        let paths = self.paths.len();
        let path_steps: usize = self.paths.iter().map(|(_, s)| s.len()).sum();
        let total_segment_bp: usize = self.node_seqs.iter().map(|s| s.len()).sum();

        // Segment length distribution.
        let mut seg_lens: Vec<usize> = self.node_seqs.iter().map(|s| s.len()).collect();
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

/// Coarse GFA topology report (V6 graph quality baseline).
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

// ── Helpers ──────────────────────────────────────────────────────

fn extract_cigar(tags: &[String]) -> Vec<CigarOp> {
    for tag in tags {
        if let Some(s) = tag.strip_prefix("cg:Z:") {
            return super::cigar::parse_cigar(s);
        }
    }
    Vec::new()
}

/// Walk a CIGAR and split at indels >= `min_var_len`, emitting segment pairs.
#[allow(clippy::too_many_arguments)]
fn split_alignment(
    tid: u32,
    t_start: i32,
    t_end: i32,
    qid: u32,
    q_start: i32,
    q_end: i32,
    q_size: i32,
    reverse: bool,
    cigar: &[CigarOp],
    min_var_len: i32,
    segments: &mut Vec<Segment>,
    links: &mut Vec<AlignmentLink>,
) {
    if cigar.is_empty() {
        // No CIGAR: treat whole record as one segment pair.
        let (qs_fwd, qe_fwd) = fwd_query_coords(q_start, q_end, q_size, reverse);
        let a = segments.len();
        segments.push(Segment {
            seq_id: tid,
            start: t_start,
            end: t_end,
        });
        let b = segments.len();
        segments.push(Segment {
            seq_id: qid,
            start: qs_fwd,
            end: qe_fwd,
        });
        links.push(AlignmentLink { a, b, reverse });
        return;
    }

    let mut ct = t_start;
    let mut cq = q_start; // in alignment-orientation coords (forward if '+', reverse if '-')
    let mut seg_start_t = ct;
    let mut seg_start_q = cq;
    let mut have_seg = false;

    let flush = |seg_start_t: i32,
                 ct: i32,
                 seg_start_q: i32,
                 cq: i32,
                 segments: &mut Vec<Segment>,
                 links: &mut Vec<AlignmentLink>,
                 have_seg: &mut bool| {
        if !*have_seg || ct <= seg_start_t {
            *have_seg = false;
            return;
        }
        let (qs_fwd, qe_fwd) = fwd_query_coords(seg_start_q, cq, q_size, reverse);
        if qs_fwd >= qe_fwd {
            *have_seg = false;
            return;
        }
        let a = segments.len();
        segments.push(Segment {
            seq_id: tid,
            start: seg_start_t,
            end: ct,
        });
        let b = segments.len();
        segments.push(Segment {
            seq_id: qid,
            start: qs_fwd,
            end: qe_fwd,
        });
        links.push(AlignmentLink { a, b, reverse });
        *have_seg = false;
    };

    for op in cigar {
        let td = op.target_delta() as i32;
        let qd = op.query_delta() as i32;
        let is_indel = (op.op() == 'I' || op.op() == 'D') && op.len() as i32 >= min_var_len;

        if is_indel {
            // Close current segment before the large indel.
            flush(
                seg_start_t,
                ct,
                seg_start_q,
                cq,
                segments,
                links,
                &mut have_seg,
            );
            ct += td;
            cq += qd;
            seg_start_t = ct;
            seg_start_q = cq;
        } else {
            // Extend current segment (small indels stay within segment as variations).
            if !have_seg {
                seg_start_t = ct;
                seg_start_q = cq;
                have_seg = true;
            }
            ct += td;
            cq += qd;
        }
    }
    // Flush trailing segment.
    flush(
        seg_start_t,
        ct,
        seg_start_q,
        cq,
        segments,
        links,
        &mut have_seg,
    );
}

/// Convert alignment-orientation query coords to forward-strand coords.
/// For '+' strand, coords are already forward. For '-' strand, flip.
fn fwd_query_coords(qs: i32, qe: i32, q_size: i32, reverse: bool) -> (i32, i32) {
    if reverse {
        (q_size - qe, q_size - qs)
    } else {
        (qs, qe)
    }
}

/// Create a novel (unaligned) node for a gap region, return its node id.
fn novel_node_for(
    node_seqs: &mut Vec<Vec<u8>>,
    sid: u32,
    start: i32,
    end: i32,
    seqs: &HashMap<String, Vec<u8>>,
    name_to_id: &HashMap<String, u32>,
) -> u32 {
    let name = id_to_name_local(name_to_id, sid).unwrap_or("?");
    let seq_bytes = seqs.get(name).map(|v| v.as_slice()).unwrap_or(&[]);
    let s = start.max(0) as usize;
    let e = (end as usize).min(seq_bytes.len());
    let bytes = if s < e {
        seq_bytes[s..e].to_vec()
    } else {
        Vec::new()
    };
    let node_id = node_seqs.len() as u32;
    node_seqs.push(bytes);
    node_id
}

fn id_to_name_local(name_to_id: &HashMap<String, u32>, id: u32) -> Option<&str> {
    name_to_id
        .iter()
        .find_map(|(name, &sid)| if sid == id { Some(name.as_str()) } else { None })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seqs_map(pairs: &[(&str, &str)]) -> HashMap<String, Vec<u8>> {
        pairs
            .iter()
            .map(|(k, v)| (k.to_string(), v.as_bytes().to_vec()))
            .collect()
    }

    #[test]
    fn test_single_alignment_no_split() {
        // 100M, no indels → one node shared by both sequences.
        let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
        let seqs = seqs_map(&[("A", &"A".repeat(100)), ("B", &"C".repeat(100))]);
        let g = PafGraph::build(paf.as_bytes(), &seqs, 100).unwrap();
        // One aligned node + possible novel trailing segments.
        // A and B should share at least one node.
        let a_nodes: Vec<u32> = g
            .paths
            .iter()
            .find(|(n, _)| n == "A")
            .unwrap()
            .1
            .iter()
            .map(|s| s.node)
            .collect();
        let b_nodes: Vec<u32> = g
            .paths
            .iter()
            .find(|(n, _)| n == "B")
            .unwrap()
            .1
            .iter()
            .map(|s| s.node)
            .collect();
        let shared = a_nodes.iter().filter(|n| b_nodes.contains(n)).count();
        assert!(shared > 0, "A and B should share a node");
    }

    #[test]
    fn test_split_at_large_indel() {
        // 50M 200I 50M: 200I >= 100 → split into two aligned nodes + one novel (insertion).
        let paf = "A\t300\t0\t100\t+\tB\t300\t0\t300\t95\t300\t255\tcg:Z:50M200I50M\n";
        let seqs = seqs_map(&[("A", &"A".repeat(300)), ("B", &"G".repeat(300))]);
        let g = PafGraph::build(paf.as_bytes(), &seqs, 100).unwrap();
        // B has an insertion of 200bp → B's path should have a novel node between aligned nodes.
        let b_path = g.paths.iter().find(|(n, _)| n == "B").unwrap();
        assert!(
            b_path.1.len() >= 3,
            "B path should have >= 3 steps (aligned, novel, aligned), got {}",
            b_path.1.len()
        );
    }

    #[test]
    fn test_small_indel_no_split() {
        // 50M 30I 50M: 30I < 100 → no split, one aligned node.
        let paf = "A\t200\t0\t130\t+\tB\t200\t0\t160\t95\t160\t255\tcg:Z:50M30I50M\n";
        let seqs = seqs_map(&[("A", &"A".repeat(200)), ("B", &"G".repeat(200))]);
        let g = PafGraph::build(paf.as_bytes(), &seqs, 100).unwrap();
        // Both A and B share exactly one aligned node for the match region.
        let a_path = g.paths.iter().find(|(n, _)| n == "A").unwrap();
        // A path: [novel 0..0? , aligned, novel trailing]. Aligned nodes should be 1.
        let aligned_in_a: Vec<u32> = a_path.1.iter().map(|s| s.node).collect();
        // The shared node (same for A and B) — check B too.
        let b_path = g.paths.iter().find(|(n, _)| n == "B").unwrap();
        let shared: Vec<u32> = aligned_in_a
            .iter()
            .filter(|n| b_path.1.iter().any(|s| &s.node == *n))
            .copied()
            .collect();
        assert_eq!(
            shared.len(),
            1,
            "exactly one shared node expected (no split), got {shared:?}"
        );
    }

    #[test]
    fn test_reverse_strand_coords_flipped() {
        // Reverse strand: query coords flipped to forward. Segments should be forward.
        let paf = "A\t100\t0\t100\t-\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
        let seqs = seqs_map(&[("A", &"A".repeat(100)), ("B", &"C".repeat(100))]);
        let g = PafGraph::build(paf.as_bytes(), &seqs, 100).unwrap();
        // Both sequences still share a node despite reverse strand.
        let a_nodes: Vec<u32> = g
            .paths
            .iter()
            .find(|(n, _)| n == "A")
            .unwrap()
            .1
            .iter()
            .map(|s| s.node)
            .collect();
        let b_nodes: Vec<u32> = g
            .paths
            .iter()
            .find(|(n, _)| n == "B")
            .unwrap()
            .1
            .iter()
            .map(|s| s.node)
            .collect();
        let shared = a_nodes.iter().filter(|n| b_nodes.contains(n)).count();
        assert!(
            shared > 0,
            "reverse-strand alignment should still produce shared node"
        );
    }

    #[test]
    fn test_gfa_output_format() {
        let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
        let seqs = seqs_map(&[("A", &"ACGT".repeat(25)), ("B", &"TGCA".repeat(25))]);
        let g = PafGraph::build(paf.as_bytes(), &seqs, 100).unwrap();
        let mut buf = Vec::new();
        g.write_gfa(&mut buf).unwrap();
        let out = String::from_utf8(buf).unwrap();
        assert!(out.starts_with("S\t1\t"), "first line should be S");
        assert!(out.contains("\nP\t"), "should contain P line");
    }
}
