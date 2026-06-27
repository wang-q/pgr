use super::index::{PafIndex, PafMetadata};
/// Disk persistence for PafIndex.
use coitrees::IntervalTree;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{self, Read, Write};
use std::path::Path;

/// File format identifier: "PGRI" = pgr index.
const MAGIC: [u8; 4] = *b"PGRI";
/// Format version (incremented on breaking changes).
const VERSION: u32 = 1;

// ── Serializable types ───────────────────────────────────────────

/// CIGAR stored as bit-packed u32 values (CigarOp.0).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlatMeta {
    pub query_id: u32,
    pub target_start: i32,
    pub target_end: i32,
    pub query_start: i32,
    pub query_end: i32,
    pub cigar: Vec<u32>,
}

/// Disk-persistable snapshot of a PafIndex.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PafIndexData {
    pub names: Vec<(String, u32)>,
    /// Per-target: (target_id, Vec<(first, last, FlatMeta)>)
    pub intervals: Vec<(u32, Vec<(i32, i32, FlatMeta)>)>,
}

// ── PafIndex serialization ───────────────────────────────────────

impl PafIndex {
    /// Save the index to a `.paf.idx` file.
    pub fn save<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let data = to_data(self);
        let encoded =
            bincode::serialize(&data).map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
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
        let _version = u32::from_le_bytes(ver_buf);
        let mut buf = Vec::new();
        f.read_to_end(&mut buf)?;
        let data: PafIndexData =
            bincode::deserialize(&buf).map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        from_data(data)
    }
}

// ── Conversion helpers ───────────────────────────────────────────

fn to_data(idx: &PafIndex) -> PafIndexData {
    let names: Vec<(String, u32)> = idx.names.iter().map(|(n, id)| (n.clone(), *id)).collect();

    let mut intervals: Vec<(u32, Vec<(i32, i32, FlatMeta)>)> = Vec::new();
    for (&tid, tree_ref) in &idx.trees {
        let mut ivs: Vec<(i32, i32, FlatMeta)> = Vec::new();
        // Iterate over the tree's intervals
        (&**tree_ref).query(
            i32::MIN,
            i32::MAX,
            |iv: &coitrees::IntervalNode<PafMetadata, u32>| {
                let m = &iv.metadata;
                ivs.push((
                    iv.first,
                    iv.last,
                    FlatMeta {
                        query_id: m.query_id,
                        target_start: m.target_start,
                        target_end: m.target_end,
                        query_start: m.query_start,
                        query_end: m.query_end,
                        cigar: m.cigar.iter().map(|op| op.0).collect(),
                    },
                ));
            },
        );
        intervals.push((tid, ivs));
    }
    PafIndexData { names, intervals }
}

fn from_data(data: PafIndexData) -> io::Result<PafIndex> {
    use coitrees::{BasicCOITree, Interval, IntervalTree};
    use indexmap::IndexMap;
    use std::collections::HashMap;
    use std::sync::Arc;

    let mut names = IndexMap::new();
    for (name, id) in &data.names {
        names.insert(name.clone(), *id);
    }

    let mut trees = HashMap::new();
    for (tid, ivs) in &data.intervals {
        let mut intervals: Vec<coitrees::Interval<PafMetadata>> = ivs
            .iter()
            .map(|(first, last, flat)| {
                let meta = PafMetadata {
                    query_id: flat.query_id,
                    target_start: flat.target_start,
                    target_end: flat.target_end,
                    query_start: flat.query_start,
                    query_end: flat.query_end,
                    cigar: flat
                        .cigar
                        .iter()
                        .map(|&v| super::cigar::CigarOp::from_raw(v))
                        .collect(),
                };
                Interval::new(*first, *last, meta)
            })
            .collect();
        intervals.sort_by(|a, b| a.first.cmp(&b.first));
        trees.insert(*tid, Arc::new(BasicCOITree::new(&intervals)));
    }

    Ok(PafIndex { names, trees })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    #[test]
    fn test_roundtrip() {
        let paf = "\
q1\t100\t0\t50\t+\tt1\t200\t0\t50\t45\t50\t255\tcg:Z:50M
q2\t300\t10\t60\t-\tt1\t200\t10\t60\t45\t50\t255\tcg:Z:50M
";
        let idx = PafIndex::build(BufReader::new(paf.as_bytes())).unwrap();
        let t1 = idx.name_to_id("t1").unwrap();

        let res_before = idx.query(t1, 0, 50);

        // Save → load
        let tmp = "/tmp/pgr_test_roundtrip.paf.idx";
        idx.save(tmp).unwrap();
        let loaded = PafIndex::load(tmp).unwrap();

        let res_after = loaded.query(t1, 0, 50);

        assert_eq!(res_before.len(), res_after.len());
        for (rb, ra) in res_before.iter().zip(res_after.iter()) {
            assert_eq!(rb.0, ra.0);
            assert_eq!(rb.1.first, ra.1.first);
            assert_eq!(rb.1.last, ra.1.last);
            assert_eq!(rb.2.first, ra.2.first);
            assert_eq!(rb.2.last, ra.2.last);
        }
    }

    #[test]
    fn test_roundtrip_with_cigar() {
        let paf = "\
A\t100\t0\t100\t+\tB\t100\t0\t100\t95\t100\t255\tcg:Z:10=5I3D80=2X
";
        let idx = PafIndex::build(BufReader::new(paf.as_bytes())).unwrap();
        let tmp = "/tmp/pgr_test_roundtrip2.paf.idx";
        idx.save(tmp).unwrap();
        let loaded = PafIndex::load(tmp).unwrap();
        let b = loaded.name_to_id("B").unwrap();
        let res = loaded.query(b, 0, 100);
        assert_eq!(res.len(), 1);
    }
}
