//! Builder for the coarse GFA graph from PAF records.

use super::{Edge, PafGraph, PathStep};
use crate::libs::paf::cigar::extract_cigar;
use crate::libs::paf::parser::parse_paf;
use std::collections::{HashMap, HashSet};
use std::io::BufRead;

use super::dsu::Dsu;
use super::segment::{novel_node_for, split_alignment, AlignmentLink, Segment};

fn flip_orient(o: char) -> char {
    if o == '+' {
        '-'
    } else {
        '+'
    }
}

/// Compute per-segment orientation relative to its node's stored sequence.
///
/// The representative segment of each DSU component (the first entry in
/// `root_info`) is oriented '+'. Each alignment link flips orientation if the
/// alignment is reverse-strand. Returns an error if a segment receives
/// conflicting orientations (e.g. an odd cycle of reverse links).
fn compute_segment_orientations(
    links: &[AlignmentLink],
    root_info: &[(usize, u32, i32, usize)],
) -> anyhow::Result<Vec<char>> {
    let max_seg_idx = links.iter().map(|l| l.a.max(l.b)).max().unwrap_or(0);
    let mut adj: Vec<Vec<(usize, bool)>> = vec![vec![]; max_seg_idx + 1];
    for link in links {
        adj[link.a].push((link.b, link.reverse));
        adj[link.b].push((link.a, link.reverse));
    }

    let mut seg_orient: Vec<Option<char>> = vec![None; adj.len()];
    for &(_, _, _, rep_idx) in root_info {
        if seg_orient[rep_idx].is_some() {
            continue;
        }
        seg_orient[rep_idx] = Some('+');
        let mut stack = vec![rep_idx];
        while let Some(curr) = stack.pop() {
            let curr_o = seg_orient[curr].expect("oriented segment");
            for &(nbr, rev) in &adj[curr] {
                let expected = if rev { flip_orient(curr_o) } else { curr_o };
                match seg_orient[nbr] {
                    Some(existing) => {
                        if existing != expected {
                            anyhow::bail!(
                                "inconsistent orientation in DSU component: segment {nbr} expected {expected} but found {existing}"
                            );
                        }
                    }
                    None => {
                        seg_orient[nbr] = Some(expected);
                        stack.push(nbr);
                    }
                }
            }
        }
    }

    Ok(seg_orient.into_iter().map(|o| o.unwrap_or('+')).collect())
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

        // Build name -> id map + reverse id -> name (reuse order from records; fall back to seqs keys).
        let mut name_to_id: HashMap<String, u32> = HashMap::new();
        let mut id_to_name: Vec<String> = Vec::new();
        let mut register = |name: &str| -> u32 {
            if let Some(&id) = name_to_id.get(name) {
                id
            } else {
                let id = id_to_name.len() as u32;
                name_to_id.insert(name.to_string(), id);
                id_to_name.push(name.to_string());
                id
            }
        };
        for r in &records {
            register(&r.target_name);
            register(&r.query_name);
        }
        if let Some(seqs_map) = seqs {
            for name in seqs_map.keys() {
                register(name);
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

        // ── Stage 2b: per-segment orientation relative to its node sequence ──
        let seg_orient = compute_segment_orientations(&links, &root_info)?;
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
            if let Some(name) = id_to_name.get(seg.seq_id as usize).map(|s| s.as_str()) {
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
        let mut seen_edges: HashSet<Edge> = HashSet::new();

        for sid in 0..num_seqs {
            let name = id_to_name
                .get(sid as usize)
                .map(|s| s.as_str())
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
                        &id_to_name,
                    );
                    steps.push(PathStep {
                        node: novel_node,
                        orient: '+',
                    });
                }
                // The aligned segment's node.
                let node = seg_node[seg_idx];
                // Orientation relative to the node's stored sequence: '+' if this
                // sequence's forward strand matches the node sequence, '-' if it
                // is the reverse complement (tracked via alignment links).
                steps.push(PathStep {
                    node,
                    orient: seg_orient[seg_idx],
                });
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
                    &id_to_name,
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
                if seen_edges.insert(e) {
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
                    &id_to_name,
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
}

#[cfg(test)]
mod tests {
    use super::{compute_segment_orientations, AlignmentLink};

    #[test]
    fn test_orientations_consistent() {
        // Two independent components:
        //   component 0: seg 0 (+) <-> seg 1 (+)  (forward link)
        //   component 1: seg 2 (+) <-> seg 3 (-)  (reverse link)
        let links = vec![
            AlignmentLink {
                a: 0,
                b: 1,
                reverse: false,
            },
            AlignmentLink {
                a: 2,
                b: 3,
                reverse: true,
            },
        ];
        let root_info = vec![(0, 0, 0, 0), (1, 1, 0, 2)];
        let orient = compute_segment_orientations(&links, &root_info).unwrap();
        assert_eq!(orient, vec!['+', '+', '+', '-']);
    }

    #[test]
    fn test_orientations_inconsistent_errors() {
        // A triangle of links where every edge is reverse-strand produces a
        // contradiction: seg 0 is forced '+' and '-' simultaneously.
        //   0 --reverse--> 1 --reverse--> 2 --reverse--> 0
        // Product of reverse flags = true * true * true = true (odd cycle).
        let links = vec![
            AlignmentLink {
                a: 0,
                b: 1,
                reverse: true,
            },
            AlignmentLink {
                a: 1,
                b: 2,
                reverse: true,
            },
            AlignmentLink {
                a: 2,
                b: 0,
                reverse: true,
            },
        ];
        let root_info = vec![(0, 0, 0, 0)];
        let err = compute_segment_orientations(&links, &root_info).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("inconsistent orientation"),
            "expected orientation error, got: {msg}"
        );
    }
}
