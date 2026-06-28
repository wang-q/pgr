/// PAF interval-tree index and query engine.
use super::cigar::{gap_compressed_identity, parse_cigar, reverse_cigar, CigarOp};
use super::parser::{parse_paf, parse_paf_line};
use super::record::PafRecord;
use coitrees::{BasicCOITree, Interval, IntervalNode, IntervalTree};
use indexmap::IndexMap;
use noodles_bgzf as bgzf;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::sync::{Arc, Mutex};

/// CIGAR storage: in-memory or lazy-loaded from a BGZF virtual position.
#[derive(Debug, Clone)]
pub enum CigarStore {
    /// CIGAR ops held in memory (in-memory build mode or after persistence).
    Owned(Vec<CigarOp>),
    /// Virtual position of the PAF line in a BGZF file (lazy fetch on demand).
    Lazy(u64),
    /// Virtual position; fetch then reverse+swap I/D (mirror entry in reverse_trees).
    LazyReversed(u64),
}

impl CigarStore {
    fn owned(cigar: Vec<CigarOp>) -> Self {
        Self::Owned(cigar)
    }
    fn lazy(vpos: u64) -> Self {
        Self::Lazy(vpos)
    }
    fn lazy_reversed(vpos: u64) -> Self {
        Self::LazyReversed(vpos)
    }
}

#[derive(Debug, Clone)]
pub struct PafMetadata {
    pub query_id: u32,
    pub target_start: i32,
    pub target_end: i32,
    pub query_start: i32,
    pub query_end: i32,
    pub strand: char,
    pub cigar: CigarStore,
}

/// Query result tuple:
/// `(query_id, query_iv, target_iv, cigar, record_target_start, record_query_start, strand)`
///
/// `record_target_start` / `record_query_start` are the original PAF record's
/// coordinates (CIGAR origin), needed for CIGAR trimming in MAF output.
/// `strand` is the original PAF record strand ('+' or '-'); needed for MAF
/// output to reverse-complement query sequences on minus-strand records.
pub type QueryResult = (
    u32,
    Interval<u32>,
    Interval<u32>,
    Vec<CigarOp>,
    i32,
    i32,
    char,
);

pub struct PafIndex {
    pub names: IndexMap<String, u32>,
    pub(crate) trees: HashMap<u32, Arc<BasicCOITree<PafMetadata, u32>>>,
    /// Mirror index: for each `+` strand record, a reversed entry is inserted
    /// into `reverse_trees[query_id]` so BFS can traverse from query → target
    /// without requiring a second PAF record with swapped roles.
    pub(crate) reverse_trees: HashMap<u32, Arc<BasicCOITree<PafMetadata, u32>>>,
    /// Lazy CIGAR source: BGZF reader + original file path (for persistence).
    pub(crate) lazy_source: Option<Mutex<bgzf::io::Reader<File>>>,
    pub(crate) lazy_source_path: Option<String>,
}

impl PafIndex {
    pub fn build<R: BufRead>(reader: R) -> std::io::Result<Self> {
        let records = parse_paf(reader)?;
        let mut names = IndexMap::new();
        let mut by_target: HashMap<u32, Vec<Interval<PafMetadata>>> = HashMap::new();
        let mut by_query: HashMap<u32, Vec<Interval<PafMetadata>>> = HashMap::new();

        for rec in &records {
            insert_record(rec, &mut names, &mut by_target, &mut by_query, None);
        }

        let trees = build_trees(by_target);
        let reverse_trees = build_trees(by_query);
        Ok(PafIndex {
            names,
            trees,
            reverse_trees,
            lazy_source: None,
            lazy_source_path: None,
        })
    }

    pub fn build_multi<R: BufRead>(readers: Vec<R>) -> std::io::Result<Self> {
        let mut names = IndexMap::new();
        let mut by_target: HashMap<u32, Vec<Interval<PafMetadata>>> = HashMap::new();
        let mut by_query: HashMap<u32, Vec<Interval<PafMetadata>>> = HashMap::new();

        for reader in readers {
            for rec in &parse_paf(reader)? {
                insert_record(rec, &mut names, &mut by_target, &mut by_query, None);
            }
        }

        let trees = build_trees(by_target);
        let reverse_trees = build_trees(by_query);
        Ok(PafIndex {
            names,
            trees,
            reverse_trees,
            lazy_source: None,
            lazy_source_path: None,
        })
    }

    /// Build an index from a file path, using lazy CIGAR loading for BGZF files.
    ///
    /// For `.gz` files that are BGZF-compressed: records the BGZF virtual position
    /// of each PAF line and stores `CigarStore::Lazy(vpos)`. CIGAR is fetched
    /// on-demand during queries, reducing memory at build time.
    ///
    /// For non-BGZF files (plain text, regular gzip): falls back to in-memory
    /// build (`CigarStore::Owned`).
    pub fn build_from_path(path: &str) -> std::io::Result<Self> {
        if path == "stdin" {
            return Self::build(crate::libs::io::reader(path));
        }

        let p = std::path::Path::new(path);
        if p.extension() == Some(std::ffi::OsStr::new("gz")) && is_bgzf(path)? {
            Self::build_lazy_bgzf(path)
        } else {
            Self::build(crate::libs::io::reader(path))
        }
    }

    fn build_lazy_bgzf(path: &str) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let mut reader = bgzf::io::Reader::new(file);

        let mut names = IndexMap::new();
        let mut by_target: HashMap<u32, Vec<Interval<PafMetadata>>> = HashMap::new();
        let mut by_query: HashMap<u32, Vec<Interval<PafMetadata>>> = HashMap::new();
        let mut line = String::new();

        loop {
            let vpos = u64::from(reader.virtual_position());
            line.clear();
            let n = reader.read_line(&mut line)?;
            if n == 0 {
                break;
            }
            let trimmed = line.trim_end_matches('\n');
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let rec = match parse_paf_line(trimmed) {
                Ok(r) => r,
                Err(e) => {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("{}: {}", e, trimmed),
                    ));
                }
            };
            insert_record(&rec, &mut names, &mut by_target, &mut by_query, Some(vpos));
        }

        let trees = build_trees(by_target);
        let reverse_trees = build_trees(by_query);
        Ok(PafIndex {
            names,
            trees,
            reverse_trees,
            lazy_source: Some(Mutex::new(bgzf::io::Reader::new(File::open(path)?))),
            lazy_source_path: Some(path.to_string()),
        })
    }

    /// Fetch CIGAR ops for a lazy record by seeking to its BGZF virtual position.
    pub fn fetch_cigar(&self, vpos: u64) -> Vec<CigarOp> {
        if let Some(ref src) = self.lazy_source {
            let mut reader = src.lock().expect("lazy_source mutex poisoned");
            if reader.seek(bgzf::VirtualPosition::from(vpos)).is_err() {
                return vec![];
            }
            let mut line = String::new();
            if reader.read_line(&mut line).is_err() {
                return vec![];
            }
            // Parse tags from the line and extract CIGAR.
            let trimmed = line.trim_end_matches('\n');
            if let Ok(rec) = parse_paf_line(trimmed) {
                return extract_cigar(&rec.tags);
            }
            vec![]
        } else {
            vec![]
        }
    }

    /// Resolve a `CigarStore` to owned CIGAR ops (clone if owned, fetch if lazy).
    fn resolve_cigar(&self, store: &CigarStore) -> Vec<CigarOp> {
        match store {
            CigarStore::Owned(c) => c.clone(),
            CigarStore::Lazy(vpos) => self.fetch_cigar(*vpos),
            CigarStore::LazyReversed(vpos) => reverse_cigar(&self.fetch_cigar(*vpos)),
        }
    }

    /// Reopen the lazy source file (used after `load` to restore lazy mode).
    pub(crate) fn reopen_lazy_source(&mut self) -> std::io::Result<()> {
        if let Some(ref path) = self.lazy_source_path {
            let file = File::open(path)?;
            self.lazy_source = Some(Mutex::new(bgzf::io::Reader::new(file)));
        }
        Ok(())
    }

    pub fn num_targets(&self) -> usize {
        self.trees.len()
    }

    /// Returns true if the index uses lazy CIGAR loading from a BGZF file.
    pub fn is_lazy(&self) -> bool {
        self.lazy_source.is_some()
    }

    pub fn name_to_id(&self, name: &str) -> Option<u32> {
        self.names.get(name).copied()
    }

    pub fn id_to_name(&self, id: u32) -> Option<&str> {
        self.names
            .get_index(id as usize)
            .map(|(name, _)| name.as_str())
    }

    pub fn query(
        &self,
        target_id: u32,
        start: i32,
        end: i32,
        min_identity: f64,
        min_output_len: i32,
    ) -> Vec<QueryResult> {
        let mut results = Vec::new();
        if let Some(tree) = self.trees.get(&target_id) {
            tree.query(start, end, |iv: &IntervalNode<PafMetadata, u32>| {
                let m = &iv.metadata;
                let cigar = self.resolve_cigar(&m.cigar);
                if gap_compressed_identity(&cigar) < min_identity {
                    return;
                }
                if let Some((qs, qe, ts, te)) = project(start, end, m, &cigar) {
                    if (qe - qs).abs() < min_output_len {
                        return;
                    }
                    results.push((
                        m.query_id,
                        Interval::new(qs, qe, m.query_id),
                        Interval::new(ts, te, target_id),
                        cigar,
                        m.target_start,
                        m.query_start,
                        m.strand,
                    ));
                }
            });
        }
        results
    }

    #[allow(clippy::too_many_arguments)]
    pub fn query_transitive_bfs(
        &self,
        target_id: u32,
        start: i32,
        end: i32,
        max_depth: u16,
        min_len: i32,
        min_dist: i32,
        min_identity: f64,
        min_output_len: i32,
        merge_distance: i32,
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
                // Query both the forward index (`trees`) and the mirror index
                // (`reverse_trees`) so BFS can traverse from any sequence in
                // both directions. Arc clones are cheap (refcount bump) and
                // release the borrow on `self` before the closure captures it.
                let fwd = self.trees.get(&tid).cloned();
                let rev = self.reverse_trees.get(&tid).cloned();
                for tree in fwd.into_iter().chain(rev) {
                    tree.query(cs, ce, |iv: &IntervalNode<PafMetadata, u32>| {
                        let m = &iv.metadata;
                        let os = cs.max(iv.first);
                        let oe = ce.min(iv.last);
                        if os >= oe {
                            return;
                        }
                        let cigar = self.resolve_cigar(&m.cigar);
                        if gap_compressed_identity(&cigar) < min_identity {
                            return;
                        }
                        if let Some((qs, qe, ts, te)) = project(os, oe, m, &cigar) {
                            if (qe - qs).abs() < min_output_len {
                                return;
                            }
                            results.push((
                                m.query_id,
                                Interval::new(qs, qe, m.query_id),
                                Interval::new(ts, te, tid),
                                cigar,
                                m.target_start,
                                m.query_start,
                                m.strand,
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
        if merge_distance > 0 {
            merge_results(&mut results, merge_distance);
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

/// Build sorted interval trees from a per-sequence interval map.
fn build_trees(
    by_seq: HashMap<u32, Vec<Interval<PafMetadata>>>,
) -> HashMap<u32, Arc<BasicCOITree<PafMetadata, u32>>> {
    let mut trees = HashMap::new();
    for (tid, mut intervals) in by_seq {
        intervals.sort_by_key(|iv| iv.first);
        trees.insert(tid, Arc::new(BasicCOITree::new(&intervals)));
    }
    trees
}

/// Insert a parsed PAF record into the name map, forward index (`by_target`),
/// and — for `+` strand records — the mirror index (`by_query`).
///
/// The mirror entry swaps query/target roles: it is stored in
/// `by_query[query_id]` with the query interval as the tree key, the
/// original target as `query_id` in metadata, and a reversed+I/D-swapped
/// CIGAR. This lets `query_transitive_bfs` traverse from any sequence in
/// both directions without requiring a second PAF record with swapped roles.
///
/// `vpos`: `Some(v)` for BGZF lazy mode (CIGAR fetched on demand),
/// `None` for in-memory mode (CIGAR stored as `Owned`).
fn insert_record(
    rec: &PafRecord,
    names: &mut IndexMap<String, u32>,
    by_target: &mut HashMap<u32, Vec<Interval<PafMetadata>>>,
    by_query: &mut HashMap<u32, Vec<Interval<PafMetadata>>>,
    vpos: Option<u64>,
) {
    let next_id = names.len() as u32;
    names.entry(rec.target_name.clone()).or_insert(next_id);
    let target_id = names[&rec.target_name];

    let next_id = names.len() as u32;
    names.entry(rec.query_name.clone()).or_insert(next_id);
    let query_id = names[&rec.query_name];

    let cigar = extract_cigar(&rec.tags);
    let (fwd_store, rev_store) = match vpos {
        Some(v) => (CigarStore::lazy(v), CigarStore::lazy_reversed(v)),
        None => (
            CigarStore::owned(cigar.clone()),
            CigarStore::owned(reverse_cigar(&cigar)),
        ),
    };

    // Forward entry: target interval → query metadata.
    // Strand is the original PAF record strand ('+' or '-'); needed for MAF
    // output to reverse-complement query sequences on minus-strand records.
    let fwd_meta = PafMetadata {
        query_id,
        target_start: rec.target_start as i32,
        target_end: rec.target_end as i32,
        query_start: rec.query_start as i32,
        query_end: rec.query_end as i32,
        strand: rec.strand,
        cigar: fwd_store,
    };
    by_target.entry(target_id).or_default().push(Interval::new(
        rec.target_start as i32,
        rec.target_end as i32,
        fwd_meta,
    ));

    // Mirror entry (reverse index): only for '+' strand records.
    // Interval is on the query coordinates; metadata.query_id is the
    // original target; query_start/end hold the original target coordinates.
    // Mirror strand is '+' because the mirror represents the query-side view
    // of an originally '+' record (query ↔ target swap preserves orientation).
    if rec.strand == '+' {
        let rev_meta = PafMetadata {
            query_id: target_id,
            target_start: rec.query_start as i32,
            target_end: rec.query_end as i32,
            query_start: rec.target_start as i32,
            query_end: rec.target_end as i32,
            strand: '+',
            cigar: rev_store,
        };
        by_query.entry(query_id).or_default().push(Interval::new(
            rec.query_start as i32,
            rec.query_end as i32,
            rev_meta,
        ));
    }
}

/// Check whether a file is BGZF-compressed by inspecting the header bytes.
fn is_bgzf(path: &str) -> std::io::Result<bool> {
    use std::io::Read;
    let mut f = File::open(path)?;
    let mut hdr = [0u8; 18];
    match f.read_exact(&mut hdr) {
        Ok(()) => {
            // BGZF: gzip magic (1f 8b 08 04), XLEN=6 at [10..12], "BC" at [12..14], SLEN=2 at [14..16]
            Ok(hdr[0] == 0x1f
                && hdr[1] == 0x8b
                && hdr[2] == 0x08
                && hdr[3] == 0x04
                && hdr[10] == 0x06
                && hdr[11] == 0x00
                && hdr[12] == b'B'
                && hdr[13] == b'C'
                && hdr[14] == 0x02
                && hdr[15] == 0x00)
        }
        Err(_) => Ok(false),
    }
}

fn project(ts: i32, te: i32, m: &PafMetadata, cigar: &[CigarOp]) -> Option<(i32, i32, i32, i32)> {
    if cigar.is_empty() {
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
    // Walk all CIGAR ops; accumulate the union of query/target intervals that
    // overlap [ts, te). Returning at the first overlap would truncate the
    // projection when the query spans indels (e.g. 4=3I3= over [0,7) must
    // project to query [0,10), not [0,4)).
    //
    // `cq` is the CIGAR query offset (0-based from the start of the aligned
    // query region). For '+' strand records this maps directly to forward
    // query coordinates (`query_start + cq`). For '-' strand records the
    // CIGAR describes RC(query) vs target, so RC offset `[rc_lo, rc_hi)` maps
    // to forward query `[query_end - rc_hi, query_end - rc_lo)`. PAF stores
    // `query_start`/`query_end` as forward-strand coordinates regardless of
    // strand, so this conversion is needed for sub-interval queries.
    let mut ct = m.target_start;
    let mut cq: i32 = 0;
    let mut q_min = i32::MAX;
    let mut q_max = i32::MIN;
    let mut t_min = i32::MAX;
    let mut t_max = i32::MIN;
    let mut found = false;
    // Convert RC offset interval [rc_lo, rc_hi) to forward query coordinates.
    let rc_to_forward = |rc_lo: i32, rc_hi: i32| -> (i32, i32) {
        if m.strand == '-' {
            (m.query_end - rc_hi, m.query_end - rc_lo)
        } else {
            (m.query_start + rc_lo, m.query_start + rc_hi)
        }
    };
    for op in cigar {
        let td = op.target_delta() as i32;
        let qd = op.query_delta() as i32;
        let ss = ct;
        let se = ct + td;
        match op.op() {
            '=' | 'X' | 'M' => {
                let os = ts.max(ss);
                let oe = te.min(se);
                if os < oe {
                    let off = os - ss;
                    let len = oe - os;
                    let (qs, qe) = rc_to_forward(cq + off, cq + off + len);
                    q_min = q_min.min(qs);
                    q_max = q_max.max(qe);
                    t_min = t_min.min(os);
                    t_max = t_max.max(oe);
                    found = true;
                }
            }
            'I' => {
                // Insertion in query at target position ct.
                // Include when ct lies within the queried target span.
                if ct >= ts && ct < te {
                    let (qs, qe) = rc_to_forward(cq, cq + qd);
                    q_min = q_min.min(qs);
                    q_max = q_max.max(qe);
                    found = true;
                }
            }
            'D' => {
                let os = ts.max(ss);
                let oe = te.min(se);
                if os < oe {
                    t_min = t_min.min(os);
                    t_max = t_max.max(oe);
                    found = true;
                }
            }
            _ => {}
        }
        ct = se;
        cq += qd;
    }
    if found {
        Some((q_min, q_max, t_min, t_max))
    } else {
        None
    }
}

fn merge_results(results: &mut Vec<QueryResult>, max_gap: i32) {
    // Group by query_id, sort by query_start, merge adjacent within max_gap
    let mut groups: HashMap<u32, Vec<(usize, i32, i32)>> = HashMap::new();
    for (i, &(qid, q_iv, _t_iv, _, _, _, _)) in results.iter().enumerate() {
        groups
            .entry(qid)
            .or_default()
            .push((i, q_iv.first, q_iv.last));
    }
    let mut to_remove = Vec::new();
    for (_qid, mut items) in groups {
        if items.len() <= 1 {
            continue;
        }
        items.sort_by_key(|&(_, s, _)| s);
        let mut prev = items[0];
        for &curr in &items[1..] {
            if (curr.1 - prev.2).abs() <= max_gap {
                // Merge: keep the one with earlier target
                results[prev.0].1.last = results[prev.0].1.last.max(curr.2);
                to_remove.push(curr.0);
            } else {
                prev = curr;
            }
        }
    }
    to_remove.sort_unstable();
    to_remove.dedup();
    for &idx in to_remove.iter().rev() {
        results.remove(idx);
    }
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

    #[test]
    fn test_merge_adjacent_intervals() {
        let mut results = vec![
            (
                0u32,
                Interval::new(0, 50, 0u32),
                Interval::new(0, 50, 1u32),
                vec![],
                0,
                0,
                '+',
            ),
            (
                0u32,
                Interval::new(55, 100, 0u32),
                Interval::new(55, 100, 1u32),
                vec![],
                55,
                55,
                '+',
            ),
        ];
        merge_results(&mut results, 10);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].1.first, 0);
        assert_eq!(results[0].1.last, 100);
    }

    #[test]
    fn test_merge_no_merge_when_far() {
        let mut results = vec![
            (
                0u32,
                Interval::new(0, 50, 0u32),
                Interval::new(0, 50, 1u32),
                vec![],
                0,
                0,
                '+',
            ),
            (
                0u32,
                Interval::new(100, 150, 0u32),
                Interval::new(100, 150, 1u32),
                vec![],
                100,
                100,
                '+',
            ),
        ];
        merge_results(&mut results, 10);
        assert_eq!(results.len(), 2);
    }
}
