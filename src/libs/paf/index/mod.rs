//! PAF interval-tree index and query engine.

mod bfs;
mod builder;
mod query;
#[cfg(test)]
mod tests;

use crate::libs::paf::cigar::{reverse_cigar, CigarOp};
use crate::libs::paf::record::PafRecord;
use coitrees::{BasicCOITree, Interval, IntervalTree};
use indexmap::IndexMap;
use noodles_bgzf as bgzf;
use std::collections::HashMap;
use std::fs::File;
use std::sync::{Arc, Mutex};

pub(crate) use query::project;

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

/// Per-record metadata stored in the interval tree, including CIGAR storage.
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

/// In-memory PAF index: name↔id mapping + per-target COI interval trees.
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
    /// Returns the number of target sequences in the index.
    pub fn num_targets(&self) -> usize {
        self.trees.len()
    }

    /// Returns true if the index uses lazy CIGAR loading from a BGZF file.
    pub fn is_lazy(&self) -> bool {
        self.lazy_source.is_some()
    }

    /// Looks up a target/query id by name.
    pub fn name_to_id(&self, name: &str) -> Option<u32> {
        self.names.get(name).copied()
    }

    /// Looks up a name by id (insertion order index).
    pub fn id_to_name(&self, id: u32) -> Option<&str> {
        self.names
            .get_index(id as usize)
            .map(|(name, _)| name.as_str())
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

    let cigar = crate::libs::paf::cigar::extract_cigar(&rec.tags)?;
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
