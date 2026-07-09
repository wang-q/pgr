/// Disk persistence for PafIndex.
use super::cigar::CigarOp;
use super::index::{CigarStore, PafIndex, PafMetadata};
use coitrees::{BasicCOITree, Interval, IntervalNode, IntervalTree};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, Write};
use std::path::Path;
use std::sync::Arc;

/// File format identifier: "PGRI" = pgr index.
const MAGIC: [u8; 4] = *b"PGRI";
/// Format version (incremented on breaking changes).
/// v3: bidirectional index — adds `reverse_intervals` and `LazyReversed`.
/// v4: per-record strand — `FlatMeta.strand` for minus-strand MAF output.
const VERSION: u32 = 4;
/// Supported format versions for deserialization.
const SUPPORTED_VERSIONS: &[u32] = &[4];

// ── Serializable types ───────────────────────────────────────────

/// CIGAR storage on disk: owned (bit-packed u32 ops), lazy (vpos), or
/// lazy-reversed (vpos + reverse/swap I/D on fetch).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum FlatCigar {
    Owned(Vec<u32>),
    Lazy(u64),
    LazyReversed(u64),
}

/// Per-record metadata stored on disk.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct FlatMeta {
    pub query_id: u32,
    pub target_start: i32,
    pub target_end: i32,
    pub query_start: i32,
    pub query_end: i32,
    pub strand: char,
    pub cigar: FlatCigar,
}

/// Flat interval list keyed by sequence id: `(seq_id, Vec<(first, last, FlatMeta)>)`.
pub type FlatTree = Vec<(u32, Vec<(i32, i32, FlatMeta)>)>;

/// Disk-persistable snapshot of a PafIndex.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct PafIndexData {
    pub names: Vec<(String, u32)>,
    /// Per-target forward intervals.
    pub intervals: FlatTree,
    /// Mirror index (reverse_trees): per-query intervals for bidirectional BFS.
    pub reverse_intervals: FlatTree,
    /// Original BGZF file path for lazy CIGAR loading (None for in-memory mode).
    pub lazy_source_path: Option<String>,
}

// ── PafIndex serialization ───────────────────────────────────────

impl PafIndex {
    /// Save the index to a `.paf.idx` file.
    pub fn save<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let data = to_data(self);
        let encoded = bincode::serialize(&data).map_err(io::Error::other)?;
        let mut f = File::create(path)?;
        f.write_all(&MAGIC)?;
        f.write_all(&VERSION.to_le_bytes())?;
        f.write_all(&encoded)?;
        Ok(())
    }

    /// Load an index from a `.paf.idx` file.
    pub fn load<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let mut f = File::open(path)?;
        let mut magic = [0u8; 4];
        f.read_exact(&mut magic)?;
        if magic != MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "not a pgr index file (bad magic)",
            ));
        }
        let mut ver_buf = [0u8; 4];
        f.read_exact(&mut ver_buf)?;
        let version = u32::from_le_bytes(ver_buf);
        if !SUPPORTED_VERSIONS.contains(&version) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "unsupported index version {version}: supported versions are {:?}; \
                     rebuild the index with `pgr paf index`",
                    SUPPORTED_VERSIONS
                ),
            ));
        }
        let mut buf = Vec::new();
        f.read_to_end(&mut buf)?;
        let data: PafIndexData = bincode::deserialize(&buf).map_err(io::Error::other)?;
        from_data(data)
    }
}

// ── Conversion helpers ───────────────────────────────────────────

/// Serialize a tree map into the flat `(seq_id, Vec<(first, last, FlatMeta)>)` form.
fn serialize_trees(trees: &HashMap<u32, Arc<BasicCOITree<PafMetadata, u32>>>) -> FlatTree {
    let mut out: FlatTree = Vec::new();
    for (&tid, tree_ref) in trees {
        let mut ivs: Vec<(i32, i32, FlatMeta)> = Vec::new();
        (**tree_ref).query(i32::MIN, i32::MAX, |iv: &IntervalNode<PafMetadata, u32>| {
            let m = &iv.metadata;
            let flat_cigar = match &m.cigar {
                CigarStore::Owned(c) => FlatCigar::Owned(c.iter().map(|op| op.0).collect()),
                CigarStore::Lazy(vpos) => FlatCigar::Lazy(*vpos),
                CigarStore::LazyReversed(vpos) => FlatCigar::LazyReversed(*vpos),
            };
            ivs.push((
                iv.first,
                iv.last,
                FlatMeta {
                    query_id: m.query_id,
                    target_start: m.target_start,
                    target_end: m.target_end,
                    query_start: m.query_start,
                    query_end: m.query_end,
                    strand: m.strand,
                    cigar: flat_cigar,
                },
            ));
        });
        out.push((tid, ivs));
    }
    out
}

/// Deserialize a flat interval list into sorted interval trees.
fn deserialize_trees(flat: &FlatTree) -> HashMap<u32, Arc<BasicCOITree<PafMetadata, u32>>> {
    use coitrees::BasicCOITree;

    let mut trees = HashMap::new();
    for (tid, ivs) in flat {
        let mut raw_intervals: Vec<Interval<PafMetadata>> = ivs
            .iter()
            .map(|(first, last, flat)| {
                let cigar = match &flat.cigar {
                    FlatCigar::Owned(vals) => {
                        CigarStore::Owned(vals.iter().map(|&v| CigarOp::from_raw(v)).collect())
                    }
                    FlatCigar::Lazy(vpos) => CigarStore::Lazy(*vpos),
                    FlatCigar::LazyReversed(vpos) => CigarStore::LazyReversed(*vpos),
                };
                let meta = PafMetadata {
                    query_id: flat.query_id,
                    target_start: flat.target_start,
                    target_end: flat.target_end,
                    query_start: flat.query_start,
                    query_end: flat.query_end,
                    strand: flat.strand,
                    cigar,
                };
                Interval::new(*first, *last, meta)
            })
            .collect();
        raw_intervals.sort_by_key(|iv| iv.first);
        trees.insert(*tid, Arc::new(BasicCOITree::new(&raw_intervals)));
    }
    trees
}

/// Export index to serializable form.
pub fn to_data(idx: &PafIndex) -> PafIndexData {
    let names: Vec<(String, u32)> = idx.names.iter().map(|(n, id)| (n.clone(), *id)).collect();
    let intervals = serialize_trees(&idx.trees);
    let reverse_intervals = serialize_trees(&idx.reverse_trees);
    PafIndexData {
        names,
        intervals,
        reverse_intervals,
        lazy_source_path: idx.lazy_source_path.clone(),
    }
}

/// Reconstruct index from serializable form.
pub fn from_data(data: PafIndexData) -> io::Result<PafIndex> {
    use indexmap::IndexMap;

    let mut names = IndexMap::new();
    for (name, id) in &data.names {
        names.insert(name.clone(), *id);
    }

    let trees = deserialize_trees(&data.intervals);
    let reverse_trees = deserialize_trees(&data.reverse_intervals);

    let mut idx = PafIndex {
        names,
        trees,
        reverse_trees,
        lazy_source: None,
        lazy_source_path: data.lazy_source_path,
    };
    // Reopen the BGZF file if a lazy source path was persisted.
    idx.reopen_lazy_source()?;
    Ok(idx)
}

// ── Tests ─────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;
    use tempfile::tempdir;

    fn build_simple() -> PafIndex {
        let paf = "\
q1\t100\t0\t50\t+\tt1\t200\t0\t50\t45\t50\t255\tcg:Z:50M
q2\t300\t10\t60\t-\tt1\t200\t10\t60\t45\t50\t255\tcg:Z:50M
";
        PafIndex::build(BufReader::new(paf.as_bytes())).unwrap()
    }

    // ── Happy path ────────────────────────────────────────────

    #[test]
    fn test_roundtrip_two_records_one_target() {
        let idx = build_simple();
        let temp = tempdir().unwrap();
        let tmp = temp.path().join("rt_2r1t.paf.idx");
        idx.save(&tmp).unwrap();
        let loaded = PafIndex::load(&tmp).unwrap();

        assert_eq!(loaded.names.len(), idx.names.len());
        assert_eq!(loaded.num_targets(), idx.num_targets());

        let t1 = loaded.name_to_id("t1").unwrap();
        let before = idx.query(t1, 0, 50, 0.0, 0);
        let after = loaded.query(t1, 0, 50, 0.0, 0);
        assert_eq!(before.len(), after.len());
        for (b, a) in before.iter().zip(after.iter()) {
            assert_eq!(b.0, a.0);
            assert_eq!(b.1.first, a.1.first);
            assert_eq!(b.1.last, a.1.last);
        }
    }

    #[test]
    fn test_roundtrip_complex_cigar() {
        let paf = "A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:10=5I3D80=2X\n";
        let idx = PafIndex::build(BufReader::new(paf.as_bytes())).unwrap();
        let temp = tempdir().unwrap();
        let tmp = temp.path().join("rt_cigar.paf.idx");
        idx.save(&tmp).unwrap();
        let loaded = PafIndex::load(&tmp).unwrap();
        let res = loaded.query(loaded.name_to_id("B").unwrap(), 0, 100, 0.0, 0);
        assert_eq!(res.len(), 1);
    }

    #[test]
    fn test_roundtrip_multi_target() {
        let paf = "\
A\t100\t0\t50\t+\tX\t200\t0\t50\t45\t50\t255\tcg:Z:50M
B\t100\t0\t50\t+\tY\t200\t0\t50\t45\t50\t255\tcg:Z:50M
C\t100\t0\t50\t+\tZ\t200\t0\t50\t45\t50\t255\tcg:Z:50M
";
        let idx = PafIndex::build(BufReader::new(paf.as_bytes())).unwrap();
        assert_eq!(idx.num_targets(), 3);

        let temp = tempdir().unwrap();
        let tmp = temp.path().join("rt_multit.paf.idx");
        idx.save(&tmp).unwrap();
        let loaded = PafIndex::load(&tmp).unwrap();

        assert_eq!(loaded.names.len(), idx.names.len());
        assert_eq!(loaded.num_targets(), 3);

        for tname in &["X", "Y", "Z"] {
            let tid = loaded.name_to_id(tname).unwrap();
            let res = loaded.query(tid, 0, 50, 0.0, 0);
            assert_eq!(res.len(), 1);
        }
    }

    #[test]
    fn test_roundtrip_transitive_bfs() {
        let paf = "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:100M
C\t100\t0\t100\t+\tA\t100\t0\t100\t90\t100\t255\tcg:Z:100M
";
        let idx = PafIndex::build(BufReader::new(paf.as_bytes())).unwrap();
        let temp = tempdir().unwrap();
        let tmp = temp.path().join("rt_bfs.paf.idx");
        idx.save(&tmp).unwrap();
        let loaded = PafIndex::load(&tmp).unwrap();

        let b = loaded.name_to_id("B").unwrap();
        let res = loaded.query_transitive_bfs(b, 0, 100, 2, 10, 10, 0.0, 0, 0, None);
        let a = loaded.name_to_id("A").unwrap();
        let c = loaded.name_to_id("C").unwrap();
        assert!(res.iter().any(|(q, _, _, _, _, _, _)| *q == a));
        assert!(res.iter().any(|(q, _, _, _, _, _, _)| *q == c));
    }

    // ── Edge cases ────────────────────────────────────────────

    #[test]
    fn test_roundtrip_empty() {
        let idx = PafIndex::build(BufReader::new("".as_bytes())).unwrap();
        assert_eq!(idx.names.len(), 0);
        assert_eq!(idx.num_targets(), 0);

        let temp = tempdir().unwrap();
        let tmp = temp.path().join("rt_empty.paf.idx");
        idx.save(&tmp).unwrap();
        let loaded = PafIndex::load(&tmp).unwrap();
        assert_eq!(loaded.names.len(), 0);
        assert_eq!(loaded.num_targets(), 0);
    }

    #[test]
    fn test_from_data_empty() {
        let data = PafIndexData {
            names: vec![],
            intervals: vec![],
            reverse_intervals: vec![],
            lazy_source_path: None,
        };
        let idx = from_data(data).unwrap();
        assert_eq!(idx.names.len(), 0);
        assert_eq!(idx.num_targets(), 0);
    }

    #[test]
    fn test_to_data_roundtrip() {
        // Verify that to_data → from_data is identity (no file I/O)
        let idx = build_simple();
        let data = to_data(&idx);
        let restored = from_data(data).unwrap();

        let t1 = restored.name_to_id("t1").unwrap();
        let res = restored.query(t1, 0, 50, 0.0, 0);
        assert_eq!(res.len(), 2);
    }

    #[test]
    fn test_to_data_preserves_names() {
        let idx = build_simple();
        let data = to_data(&idx);
        assert_eq!(data.names.len(), 3); // q1, q2, t1
        let names: Vec<&str> = data.names.iter().map(|(n, _)| n.as_str()).collect();
        assert!(names.contains(&"q1"));
        assert!(names.contains(&"q2"));
        assert!(names.contains(&"t1"));
    }

    // ── Error handling ────────────────────────────────────────

    #[test]
    fn test_load_bad_magic() {
        let temp = tempdir().unwrap();
        let tmp = temp.path().join("bad_magic.paf.idx");
        std::fs::write(&tmp, "XXXXsome garbage").unwrap();
        let err = PafIndex::load(&tmp).err().unwrap();
        assert!(err.to_string().contains("bad magic"));
    }

    #[test]
    fn test_load_truncated_header() {
        let temp = tempdir().unwrap();
        let tmp = temp.path().join("trunc_hdr.paf.idx");
        std::fs::write(&tmp, "PG").unwrap(); // only 2 bytes, can't read 4-byte header
        assert!(PafIndex::load(&tmp).is_err());
    }

    #[test]
    fn test_load_truncated_body() {
        let temp = tempdir().unwrap();
        let tmp = temp.path().join("trunc_body.paf.idx");
        let mut f = std::fs::File::create(&tmp).unwrap();
        f.write_all(&MAGIC).unwrap();
        f.write_all(&VERSION.to_le_bytes()).unwrap();
        f.write_all(b"incomplete bincode").unwrap();
        // Should fail during bincode deserialization
        assert!(PafIndex::load(&tmp).is_err());
    }

    #[test]
    fn test_save_load_idempotent() {
        // save → load → save → load: second round should match first
        let idx = build_simple();
        let temp = tempdir().unwrap();
        let tmp1 = temp.path().join("idem_1.paf.idx");
        let tmp2 = temp.path().join("idem_2.paf.idx");

        idx.save(&tmp1).unwrap();
        let loaded1 = PafIndex::load(&tmp1).unwrap();
        loaded1.save(&tmp2).unwrap();
        let loaded2 = PafIndex::load(&tmp2).unwrap();

        let t1 = loaded2.name_to_id("t1").unwrap();
        assert_eq!(loaded2.query(t1, 0, 50, 0.0, 0).len(), 2);
    }
}
