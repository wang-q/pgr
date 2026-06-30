/// PAF interval-tree index and query engine.
use super::cigar::{extract_cigar, gap_compressed_identity, reverse_cigar, CigarOp};
use super::parser::{parse_paf, parse_paf_line};
use super::record::PafRecord;
use coitrees::{BasicCOITree, Interval, IntervalNode, IntervalTree};
use indexmap::IndexMap;
use noodles_bgzf as bgzf;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::sync::{Arc, Mutex};

mod bfs;
#[cfg(test)]
mod tests;

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
            insert_record(rec, &mut names, &mut by_target, &mut by_query, None)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()))?;
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
                insert_record(rec, &mut names, &mut by_target, &mut by_query, None).map_err(
                    |e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()),
                )?;
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
    pub fn build_from_path(path: &str) -> anyhow::Result<Self> {
        if path == "stdin" {
            return Ok(Self::build(crate::libs::io::reader(path)?)?);
        }

        let p = std::path::Path::new(path);
        if p.extension() == Some(std::ffi::OsStr::new("gz")) && crate::is_bgzf(path) {
            Ok(Self::build_lazy_bgzf(path)?)
        } else {
            Ok(Self::build(crate::libs::io::reader(path)?)?)
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
            insert_record(&rec, &mut names, &mut by_target, &mut by_query, Some(vpos))
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()))?;
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
                return extract_cigar(&rec.tags).unwrap_or_default();
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
) -> anyhow::Result<()> {
    let next_id = names.len() as u32;
    names.entry(rec.target_name.clone()).or_insert(next_id);
    let target_id = names[&rec.target_name];

    let next_id = names.len() as u32;
    names.entry(rec.query_name.clone()).or_insert(next_id);
    let query_id = names[&rec.query_name];

    let cigar = extract_cigar(&rec.tags)?;
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

    Ok(())
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
