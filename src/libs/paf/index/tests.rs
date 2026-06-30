use super::*;
use crate::libs::paf::cigar::extract_cigar;
use std::io::BufReader;

fn paf_data() -> &'static str {
    "\
q1\t100\t0\t50\t+\tt1\t200\t0\t50\t45\t50\t255\tcg:Z:50M\tgi:f:0.9
q2\t300\t10\t60\t-\tt1\t200\t10\t60\t45\t50\t255\tcg:Z:50M\tgi:f:0.9
q3\t400\t0\t40\t+\tt2\t500\t0\t40\t38\t40\t255\tcg:Z:40M
"
}

#[test]
fn test_build() {
    let idx = PafIndex::build(BufReader::new(paf_data().as_bytes())).unwrap();
    assert_eq!(idx.names.len(), 5);
    assert_eq!(idx.num_targets(), 2);
}

#[test]
fn test_query() {
    let idx = PafIndex::build(BufReader::new(paf_data().as_bytes())).unwrap();
    let t1 = idx.name_to_id("t1").unwrap();
    let res = idx.query(t1, 0, 50, 0.0, 0);
    assert_eq!(res.len(), 2, "expected 2 overlapping records for t1:[0,50)");
    let qids: Vec<u32> = res.iter().map(|(q, _, _, _, _, _, _)| *q).collect();
    assert!(
        qids.contains(&idx.name_to_id("q1").unwrap()),
        "q1 not found"
    );
    assert!(
        qids.contains(&idx.name_to_id("q2").unwrap()),
        "q2 not found"
    );
    assert_eq!(res[0].0, idx.name_to_id("q1").unwrap());
}

#[test]
fn test_query_no_overlap() {
    let idx = PafIndex::build(BufReader::new(paf_data().as_bytes())).unwrap();
    let t1 = idx.name_to_id("t1").unwrap();
    assert!(idx.query(t1, 100, 150, 0.0, 0).is_empty());
}

#[test]
fn test_bfs_two_hop() {
    let paf = "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t100\t+\tA\t100\t0\t100\t90\t100\t255\tcg:Z:100M
";
    let idx = PafIndex::build(BufReader::new(paf.as_bytes())).unwrap();
    let b = idx.name_to_id("B").unwrap();
    let res = idx.query_transitive_bfs(b, 0, 100, 2, 10, 10, 0.0, 0, 0);
    let a = idx.name_to_id("A").unwrap();
    let c = idx.name_to_id("C").unwrap();
    assert!(
        res.iter().any(|(q, _, _, _, _, _, _)| *q == a),
        "A not found"
    );
    assert!(
        res.iter().any(|(q, _, _, _, _, _, _)| *q == c),
        "C not found"
    );
}

#[test]
fn test_extract_cigar() {
    let c = extract_cigar(&["cg:Z:10=5I3D".into(), "gi:f:0.9".into()]).unwrap();
    assert_eq!(c.len(), 3);
}

#[test]
fn test_extract_cigar_empty() {
    assert!(extract_cigar(&["gi:f:0.9".into()]).unwrap().is_empty());
}

#[test]
fn test_build_multi_merges_targets() {
    let paf1 = "A\t100\t0\t50\t+\tX\t200\t0\t50\t45\t50\t255\tcg:Z:50M\n";
    let paf2 = "B\t100\t0\t50\t+\tX\t200\t50\t100\t45\t50\t255\tcg:Z:50M\n";
    let idx = PafIndex::build_multi(vec![
        BufReader::new(paf1.as_bytes()),
        BufReader::new(paf2.as_bytes()),
    ])
    .unwrap();
    // X is shared target across both files → 1 target
    assert_eq!(idx.num_targets(), 1);
    assert_eq!(idx.names.len(), 3); // A, B, X
    let x = idx.name_to_id("X").unwrap();
    let res = idx.query(x, 0, 100, 0.0, 0);
    assert_eq!(res.len(), 2);
}

// ── project() edge cases ──────────────────────────────────

#[test]
fn test_project_empty_cigar_outside() {
    let m = PafMetadata {
        query_id: 0,
        target_start: 0,
        target_end: 50,
        query_start: 0,
        query_end: 50,
        strand: '+',
        cigar: CigarStore::owned(vec![]),
    };
    assert!(project(100, 200, &m, &[]).is_none());
}

#[test]
fn test_project_cigar_no_overlap() {
    let cigar = vec![CigarOp::new(50, 'M')];
    let m = PafMetadata {
        query_id: 0,
        target_start: 0,
        target_end: 50,
        query_start: 0,
        query_end: 50,
        strand: '+',
        cigar: CigarStore::owned(cigar.clone()),
    };
    assert!(project(100, 200, &m, &cigar).is_none());
}

#[test]
fn test_project_cigar_with_insertion_offset() {
    // CIGAR: 10M5I10M. Query [11,16) on target lands in the trailing M segment,
    // but query coordinates are shifted by the 5-base insertion.
    let cigar = vec![
        CigarOp::new(10, 'M'),
        CigarOp::new(5, 'I'),
        CigarOp::new(10, 'M'),
    ];
    let m = PafMetadata {
        query_id: 0,
        target_start: 0,
        target_end: 25,
        query_start: 0,
        query_end: 25,
        strand: '+',
        cigar: CigarStore::owned(cigar.clone()),
    };
    let (qs, qe, ts, te) = project(11, 16, &m, &cigar).unwrap();
    assert_eq!((qs, qe, ts, te), (16, 21, 11, 16));
}

// ── project() on '-' strand sub-intervals ────────────────
//
// PAF '-' strand: query_start/query_end are forward-strand coordinates,
// but CIGAR describes RC(query) vs target. RC offset [rc_lo, rc_hi) maps
// to forward [query_end - rc_hi, query_end - rc_lo). Full-overlap queries
// return the entire forward region (strand-agnostic), but sub-intervals
// must reverse the offset mapping.

fn minus_metadata(qs: i32, qe: i32, ts: i32, te: i32) -> PafMetadata {
    PafMetadata {
        query_id: 0,
        target_start: ts,
        target_end: te,
        query_start: qs,
        query_end: qe,
        strand: '-',
        cigar: CigarStore::owned(vec![]),
    }
}

#[test]
fn test_project_minus_strand_full_overlap() {
    // 10= over forward query [0,10); query the full target [0,10).
    // RC offset [0,10) → forward [10-10, 10-0) = [0,10) — same as full
    // record, so '+' and '-' strand agree here.
    let cigar = vec![CigarOp::new(10, '=')];
    let m = minus_metadata(0, 10, 0, 10);
    let (qs, qe, ts, te) = project(0, 10, &m, &cigar).unwrap();
    assert_eq!((qs, qe, ts, te), (0, 10, 0, 10));
}

#[test]
fn test_project_minus_strand_subinterval_first_half() {
    // 10= over forward query [0,10). Query target [0,5) overlaps the
    // first 5 CIGAR query bases = RC(query)[0..5] = complement of
    // query[5..10] reversed → forward [5,10).
    let cigar = vec![CigarOp::new(10, '=')];
    let m = minus_metadata(0, 10, 0, 10);
    let (qs, qe, ts, te) = project(0, 5, &m, &cigar).unwrap();
    assert_eq!((qs, qe, ts, te), (5, 10, 0, 5));
}

#[test]
fn test_project_minus_strand_subinterval_second_half() {
    // 10= over forward query [0,10). Query target [5,10) overlaps the
    // last 5 CIGAR query bases = RC(query)[5..10] = complement of
    // query[0..5] reversed → forward [0,5).
    let cigar = vec![CigarOp::new(10, '=')];
    let m = minus_metadata(0, 10, 0, 10);
    let (qs, qe, ts, te) = project(5, 10, &m, &cigar).unwrap();
    assert_eq!((qs, qe, ts, te), (0, 5, 5, 10));
}

#[test]
fn test_project_minus_strand_with_query_offset() {
    // 10= over forward query [100,110). Full overlap → forward [100,110).
    let cigar = vec![CigarOp::new(10, '=')];
    let m = minus_metadata(100, 110, 0, 10);
    let (qs, qe, _ts, _te) = project(0, 10, &m, &cigar).unwrap();
    assert_eq!((qs, qe), (100, 110));
    // Sub-interval target [0,5) → forward [105,110).
    let (qs, qe, _, _) = project(0, 5, &m, &cigar).unwrap();
    assert_eq!((qs, qe), (105, 110));
    // Sub-interval target [5,10) → forward [100,105).
    let (qs, qe, _, _) = project(5, 10, &m, &cigar).unwrap();
    assert_eq!((qs, qe), (100, 105));
}

#[test]
fn test_project_minus_strand_with_insertion() {
    // CIGAR: 5=3I2= over forward query [0,10). Target span = 7.
    // RC offset walk: op1 5= covers RC[0,5); op2 3I at RC[5,8); op3 2=
    // covers RC[8,10).
    // Query target [0,5) hits op1 only (op2 sits at target pos 5, outside
    // the half-open [0,5)) → forward [10-5, 10-0) = [5,10).
    // Query target [5,7) hits op2 (insertion at target pos 5) AND op3:
    //   op2 RC[5,8) → forward [2,5); op3 RC[8,10) → forward [0,2)
    //   union = forward [0,5).
    // Query target [0,7) hits all three ops → forward [0,10).
    let cigar = vec![
        CigarOp::new(5, '='),
        CigarOp::new(3, 'I'),
        CigarOp::new(2, '='),
    ];
    let m = minus_metadata(0, 10, 0, 7);
    let (qs, qe, _, _) = project(0, 5, &m, &cigar).unwrap();
    assert_eq!((qs, qe), (5, 10));
    let (qs, qe, _, _) = project(5, 7, &m, &cigar).unwrap();
    assert_eq!((qs, qe), (0, 5));
    let (qs, qe, _, _) = project(0, 7, &m, &cigar).unwrap();
    assert_eq!((qs, qe), (0, 10));
}

#[test]
fn test_query_min_identity_filters() {
    let idx = PafIndex::build(BufReader::new(paf_data().as_bytes())).unwrap();
    let t1 = idx.name_to_id("t1").unwrap();
    let res = idx.query(t1, 0, 50, 0.95, 0);
    assert_eq!(res.len(), 2);
    let res = idx.query(t1, 0, 50, 1.01, 0);
    assert_eq!(res.len(), 0);
}
