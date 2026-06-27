/// PAF interval-tree index and query engine.
use super::cigar::{parse_cigar, CigarOp};
use super::parser::parse_paf;
use coitrees::{BasicCOITree, Interval, IntervalNode, IntervalTree};
use indexmap::IndexMap;
use std::collections::HashMap;
use std::io::BufRead;
use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct PafMetadata {
    pub query_id: u32,
    pub target_start: i32,
    pub target_end: i32,
    pub query_start: i32,
    pub query_end: i32,
    pub cigar: Vec<CigarOp>,
}

pub type QueryResult = (u32, Interval<u32>, Interval<u32>);

pub struct PafIndex {
    pub names: IndexMap<String, u32>,
    pub(crate) trees: HashMap<u32, Arc<BasicCOITree<PafMetadata, u32>>>,
}

impl PafIndex {
    pub fn build<R: BufRead>(reader: R) -> std::io::Result<Self> {
        let records = parse_paf(reader)?;
        let mut names = IndexMap::new();
        let mut by_target: HashMap<u32, Vec<Interval<PafMetadata>>> = HashMap::new();

        for rec in &records {
            let next_id = names.len() as u32;
            names.entry(rec.target_name.clone()).or_insert(next_id);
            let target_id = names[&rec.target_name];

            let next_id = names.len() as u32;
            names.entry(rec.query_name.clone()).or_insert(next_id);
            let query_id = names[&rec.query_name];

            let cigar = extract_cigar(&rec.tags);
            let meta = PafMetadata {
                query_id,
                target_start: rec.target_start as i32,
                target_end: rec.target_end as i32,
                query_start: rec.query_start as i32,
                query_end: rec.query_end as i32,
                cigar,
            };
            by_target.entry(target_id).or_default().push(Interval::new(
                rec.target_start as i32,
                rec.target_end as i32,
                meta,
            ));
        }

        let mut trees = HashMap::new();
        for (tid, mut intervals) in by_target {
            intervals.sort_by(|a, b| a.first.cmp(&b.first));
            trees.insert(tid, Arc::new(BasicCOITree::new(&intervals)));
        }
        Ok(PafIndex { names, trees })
    }

    pub fn num_targets(&self) -> usize {
        self.trees.len()
    }

    pub fn name_to_id(&self, name: &str) -> Option<u32> {
        self.names.get(name).copied()
    }

    pub fn id_to_name(&self, id: u32) -> Option<&str> {
        self.names
            .get_index(id as usize)
            .map(|(name, _)| name.as_str())
    }

    pub fn query(&self, target_id: u32, start: i32, end: i32) -> Vec<QueryResult> {
        let mut results = Vec::new();
        if let Some(tree) = self.trees.get(&target_id) {
            tree.query(start, end, |iv: &IntervalNode<PafMetadata, u32>| {
                let m = &iv.metadata;
                if let Some((qs, qe, ts, te)) = project(start, end, m) {
                    results.push((
                        m.query_id,
                        Interval::new(qs, qe, m.query_id),
                        Interval::new(ts, te, target_id),
                    ));
                }
            });
        }
        results
    }

    pub fn query_transitive_bfs(
        &self,
        target_id: u32,
        start: i32,
        end: i32,
        max_depth: u16,
        min_len: i32,
        min_dist: i32,
    ) -> Vec<QueryResult> {
        let mut results = Vec::new();
        let mut visited: HashMap<u32, SortedRanges> = HashMap::new();

        let init = visited
            .entry(target_id)
            .or_insert_with(|| SortedRanges::new(i32::MAX, min_dist))
            .insert(start, end);
        let mut current: Vec<(u32, i32, i32)> = init
            .into_iter()
            .filter(|&(s, e)| (e - s).abs() >= min_len)
            .map(|(s, e)| (target_id, s, e))
            .collect();

        let mut depth = 0u16;
        while !current.is_empty() && (max_depth == 0 || depth < max_depth) {
            let mut next = vec![];
            for &(tid, cs, ce) in &current {
                if let Some(tree) = self.trees.get(&tid) {
                    tree.query(cs, ce, |iv: &IntervalNode<PafMetadata, u32>| {
                        let m = &iv.metadata;
                        let os = cs.max(iv.first);
                        let oe = ce.min(iv.last);
                        if os >= oe {
                            return;
                        }
                        if let Some((qs, qe, ts, te)) = project(os, oe, m) {
                            results.push((
                                m.query_id,
                                Interval::new(qs, qe, m.query_id),
                                Interval::new(ts, te, tid),
                            ));
                            if m.query_id != tid {
                                let sr = visited
                                    .entry(m.query_id)
                                    .or_insert_with(|| SortedRanges::new(i32::MAX, min_dist));
                                for (ns, ne) in sr.insert(qs, qe) {
                                    if (ne - ns).abs() >= min_len {
                                        next.push((m.query_id, ns, ne));
                                    }
                                }
                            }
                        }
                    });
                }
            }
            current = next;
            depth += 1;
        }
        results
    }
}

fn extract_cigar(tags: &[String]) -> Vec<CigarOp> {
    for tag in tags {
        if let Some(s) = tag.strip_prefix("cg:Z:") {
            return parse_cigar(s);
        }
    }
    vec![]
}

fn project(ts: i32, te: i32, m: &PafMetadata) -> Option<(i32, i32, i32, i32)> {
    if m.cigar.is_empty() {
        let off = (ts - m.target_start).max(0);
        let len = (te - ts).min((m.target_end - m.target_start) - off);
        if len <= 0 {
            return None;
        }
        return Some((
            m.query_start + off,
            m.query_start + off + len,
            ts.max(m.target_start),
            (ts + len).min(m.target_end),
        ));
    }
    let mut ct = m.target_start;
    let mut cq = m.query_start;
    for op in &m.cigar {
        let td = op.target_delta() as i32;
        let qd = op.query_delta() as i32;
        let ss = ct;
        let se = ct + td;
        let os = ts.max(ss);
        let oe = te.min(se);
        if os < oe {
            let off = os - ss;
            let len = oe - os;
            let qoff = if op.op() == 'I' { 0 } else { off };
            return Some((cq + qoff, cq + qoff + len, os, oe));
        }
        ct = se;
        cq += qd;
    }
    None
}

struct SortedRanges {
    ranges: Vec<(i32, i32)>,
    seq_len: i32,
    min_dist: i32,
}

impl SortedRanges {
    fn new(seq_len: i32, min_dist: i32) -> Self {
        Self {
            ranges: vec![],
            seq_len,
            min_dist,
        }
    }

    fn insert(&mut self, start: i32, end: i32) -> Vec<(i32, i32)> {
        let (mut s, mut e) = if start <= end {
            (start, end)
        } else {
            (end, start)
        };
        s = s.max(0);
        e = e.min(self.seq_len);
        if s >= e {
            return vec![];
        }
        // Find insertion point
        let mut i = match self.ranges.binary_search_by_key(&s, |&(a, _)| a) {
            Ok(p) => p,
            Err(p) => p,
        };
        // Step back if previous range overlaps
        if i > 0 && self.ranges[i - 1].1 > s {
            i -= 1;
        }
        // Expand to nearby ranges
        if i > 0 && (s - self.ranges[i - 1].1).abs() <= self.min_dist {
            s = self.ranges[i - 1].1;
            i -= 1;
        }
        if i < self.ranges.len() && (self.ranges[i].0 - e).abs() <= self.min_dist {
            e = self.ranges[i].0;
        }
        // Collect non-covered portions
        let mut out = vec![];
        let mut cur = s;
        let mut j = i;
        while j < self.ranges.len() && cur < e {
            let (rs, re) = self.ranges[j];
            if rs > e {
                break;
            }
            if cur < rs {
                out.push((cur, rs.min(e)));
            }
            cur = cur.max(re);
            j += 1;
        }
        if cur < e {
            out.push((cur, e));
        }
        // Merge
        let first = match self.ranges.binary_search_by_key(&s, &|&(a, _)| a) {
            Ok(p) => p,
            Err(p) => {
                if p > 0 && self.ranges[p - 1].1 >= s {
                    p - 1
                } else {
                    p
                }
            }
        };
        let last = self.ranges[first..]
            .iter()
            .position(|&(a, _)| a > e)
            .map(|p| first + p)
            .unwrap_or(self.ranges.len());
        let old_end = self.ranges.get(last.saturating_sub(1)).map(|r| r.1);
        let merged = (
            s.min(self.ranges.get(first).map(|r| r.0).unwrap_or(s)),
            e.max(old_end.unwrap_or(e)),
        );
        self.ranges.drain(first..last);
        self.ranges.insert(first, merged);
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
        let res = idx.query(t1, 0, 50);
        assert_eq!(res.len(), 2, "expected 2 overlapping records for t1:[0,50)");
        let qids: Vec<u32> = res.iter().map(|(q, _, _)| *q).collect();
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
        assert!(idx.query(t1, 100, 150).is_empty());
    }

    #[test]
    fn test_bfs_two_hop() {
        let paf = "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t100\t+\tA\t100\t0\t100\t90\t100\t255\tcg:Z:100M
";
        let idx = PafIndex::build(BufReader::new(paf.as_bytes())).unwrap();
        let b = idx.name_to_id("B").unwrap();
        let res = idx.query_transitive_bfs(b, 0, 100, 2, 10, 10);
        let a = idx.name_to_id("A").unwrap();
        let c = idx.name_to_id("C").unwrap();
        assert!(res.iter().any(|(q, _, _)| *q == a), "A not found");
        assert!(res.iter().any(|(q, _, _)| *q == c), "C not found");
    }

    #[test]
    fn test_sorted_ranges_disjoint() {
        let mut sr = SortedRanges::new(1000, 0);
        assert_eq!(sr.insert(100, 200), vec![(100, 200)]);
        assert_eq!(sr.insert(300, 400), vec![(300, 400)]);
    }

    #[test]
    fn test_sorted_ranges_contained() {
        let mut sr = SortedRanges::new(1000, 0);
        sr.insert(100, 200);
        assert!(sr.insert(120, 180).is_empty());
    }

    #[test]
    fn test_sorted_ranges_adjacent() {
        let mut sr = SortedRanges::new(1000, 10);
        sr.insert(100, 200);
        sr.insert(210, 300);
        assert_eq!(sr.ranges.len(), 1);
        assert_eq!(sr.ranges[0], (100, 300));
    }

    #[test]
    fn test_extract_cigar() {
        let c = extract_cigar(&["cg:Z:10=5I3D".into(), "gi:f:0.9".into()]);
        assert_eq!(c.len(), 3);
    }

    #[test]
    fn test_extract_cigar_empty() {
        assert!(extract_cigar(&["gi:f:0.9".into()]).is_empty());
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
            cigar: vec![],
        };
        assert!(project(100, 200, &m).is_none());
    }

    #[test]
    fn test_project_cigar_no_overlap() {
        let m = PafMetadata {
            query_id: 0,
            target_start: 0,
            target_end: 50,
            query_start: 0,
            query_end: 50,
            cigar: vec![CigarOp::new(50, 'M')],
        };
        assert!(project(100, 200, &m).is_none());
    }

    #[test]
    fn test_project_cigar_with_insertion_offset() {
        // CIGAR: 10M5I10M. Query [11,16) on target lands in the trailing M segment,
        // but query coordinates are shifted by the 5-base insertion.
        let m = PafMetadata {
            query_id: 0,
            target_start: 0,
            target_end: 25,
            query_start: 0,
            query_end: 25,
            cigar: vec![
                CigarOp::new(10, 'M'),
                CigarOp::new(5, 'I'),
                CigarOp::new(10, 'M'),
            ],
        };
        let (qs, qe, ts, te) = project(11, 16, &m).unwrap();
        assert_eq!((qs, qe, ts, te), (16, 21, 11, 16));
    }
}
