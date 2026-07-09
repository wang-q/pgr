use super::*;
use std::collections::HashMap;

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
    let g = PafGraph::build(paf.as_bytes(), Some(&seqs), 100).unwrap();
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
    let g = PafGraph::build(paf.as_bytes(), Some(&seqs), 100).unwrap();
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
    let g = PafGraph::build(paf.as_bytes(), Some(&seqs), 100).unwrap();
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
    let g = PafGraph::build(paf.as_bytes(), Some(&seqs), 100).unwrap();
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
fn test_reverse_strand_path_orientation() {
    // A = ACGTACGTAC, B = GTACGTACGT (reverse complement of A).
    // PAF order: query first, target second. B is the query aligned in reverse
    // to target A; they should share one node, but B's path step must traverse
    // it in the '-' orientation.
    let paf = "B\t10\t0\t10\t-\tA\t10\t0\t10\t10\t10\t255\tcg:Z:10M\n";
    let seqs = seqs_map(&[("A", "ACGTACGTAC"), ("B", "GTACGTACGT")]);
    let g = PafGraph::build(paf.as_bytes(), Some(&seqs), 100).unwrap();
    let a_path = g.paths.iter().find(|(n, _)| n == "A").unwrap();
    let b_path = g.paths.iter().find(|(n, _)| n == "B").unwrap();
    assert_eq!(a_path.1.len(), 1, "A should have one aligned step");
    assert_eq!(b_path.1.len(), 1, "B should have one aligned step");
    assert_eq!(
        a_path.1[0].node, b_path.1[0].node,
        "A and B should share a node"
    );
    assert_eq!(a_path.1[0].orient, '+', "A traverses the node forward");
    assert_eq!(b_path.1[0].orient, '-', "B traverses the node reverse");
    assert_eq!(
        String::from_utf8_lossy(&g.node_seqs[a_path.1[0].node as usize]),
        "ACGTACGTAC",
        "node sequence should match A's forward strand"
    );
}

#[test]
fn test_gfa_output_format() {
    let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M\n";
    let seqs = seqs_map(&[("A", &"ACGT".repeat(25)), ("B", &"TGCA".repeat(25))]);
    let g = PafGraph::build(paf.as_bytes(), Some(&seqs), 100).unwrap();
    let mut buf = Vec::new();
    g.write_gfa(&mut buf).unwrap();
    let out = String::from_utf8(buf).unwrap();
    assert!(out.starts_with("S\t1\t"), "first line should be S");
    assert!(out.contains("\nP\t"), "should contain P line");
}
