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
