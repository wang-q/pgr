//! Builder for the coarse GFA graph from PAF records.

use super::{Edge, PafGraph, PathStep};
use crate::libs::paf::cigar::extract_cigar;
use crate::libs::paf::parser::parse_paf;
use std::collections::HashMap;
use std::io::BufRead;

use super::dsu::Dsu;
use super::segment::{id_to_name_local, novel_node_for, split_alignment, AlignmentLink, Segment};

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
}
