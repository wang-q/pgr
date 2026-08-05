//! Project a `.pgi` index onto a hypervector for cheap distance comparisons.

use super::{count_unique, PgiQuery};
use anyhow::Context;
use std::io::{Read, Write};

/// Hypervector file magic.
pub const HV_MAGIC: &[u8; 4] = b"PGV1";
/// Hypervector file version (v2: sparse projection, stores s and k-mer count).
pub const HV_VERSION: u32 = 2;

/// Fold a u128 k-mer key into a u64 seed for HV projection.
fn key_to_seed(kmer: u128) -> u64 {
    (kmer as u64) ^ ((kmer >> 64) as u64)
}

/// Project the index's unique k-mer keys onto a sparse `dim`-dimension
/// hypervector, each key updating `sparse` random dimensions with ±1.
pub fn index_to_hv(idx: &impl PgiQuery, dim: usize, sparse: usize) -> Vec<i32> {
    let (i0, i1) = idx.entry_range(0, u128::MAX);
    let n = count_unique(idx) as usize;
    let mut seeds = Vec::with_capacity(n);
    let mut i = i0;
    while i < i1 {
        seeds.push(key_to_seed(idx.entry_kmer(i)));
        i = idx.entry_next(i);
    }
    crate::libs::hv::hash_hv_sparse(&seeds, dim, sparse)
}

/// Serialize a hypervector to the `.hv` file format.
pub fn write_hv<W: Write>(
    w: &mut W,
    name: &str,
    k: usize,
    dim: usize,
    sparse: usize,
    n_kmer: usize,
    hv: &[i32],
) -> anyhow::Result<()> {
    anyhow::ensure!(hv.len() == dim, "hv length mismatch");
    w.write_all(HV_MAGIC)?;
    w.write_all(&HV_VERSION.to_le_bytes())?;
    w.write_all(&(k as u32).to_le_bytes())?;
    w.write_all(&(dim as u32).to_le_bytes())?;
    w.write_all(&(sparse as u32).to_le_bytes())?;
    w.write_all(&(n_kmer as u64).to_le_bytes())?;
    let nb = name.len() as u32;
    w.write_all(&nb.to_le_bytes())?;
    w.write_all(name.as_bytes())?;
    for v in hv {
        w.write_all(&v.to_le_bytes())?;
    }
    Ok(())
}

/// Deserialize a hypervector from the `.hv` file format.
pub struct HvFile {
    pub name: String,
    pub k: usize,
    pub dim: usize,
    /// Dimensions updated per k-mer by the sparse projection.
    pub sparse: usize,
    /// Number of unique k-mers projected (set cardinality).
    pub n_kmer: usize,
    pub hv: Vec<i32>,
}

/// Load a `.hv` file.
pub fn read_hv<R: Read>(r: &mut R) -> anyhow::Result<HvFile> {
    let mut magic = [0u8; 4];
    r.read_exact(&mut magic).context("reading hv magic")?;
    if &magic != HV_MAGIC {
        anyhow::bail!("not a pgr hv file (bad magic)");
    }
    let version = read_u32(r)?;
    if version != HV_VERSION {
        anyhow::bail!("unsupported hv version {version}");
    }
    let k = read_u32(r)? as usize;
    let dim = read_u32(r)? as usize;
    let sparse = read_u32(r)? as usize;
    let n_kmer = read_u64(r)? as usize;
    let nb = read_u32(r)? as usize;
    let mut name = Vec::new();
    name.try_reserve_exact(nb).context("hv name too large")?;
    name.resize(nb, 0);
    r.read_exact(&mut name)?;
    let name = String::from_utf8(name).context("hv name utf8")?;
    let mut hv = Vec::new();
    hv.try_reserve_exact(dim)
        .context("hv dimension too large")?;
    // Read element by element: a crafted `dim` is rejected as soon as the
    // input runs out, without first zero-filling a huge allocation.
    let mut buf = [0u8; 4];
    for _ in 0..dim {
        r.read_exact(&mut buf)?;
        hv.push(i32::from_le_bytes(buf));
    }
    Ok(HvFile {
        name,
        k,
        dim,
        sparse,
        n_kmer,
        hv,
    })
}

fn read_u32<R: Read>(r: &mut R) -> anyhow::Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

fn read_u64<R: Read>(r: &mut R) -> anyhow::Result<u64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(u64::from_le_bytes(b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::pgi::build::build_from_seqs;

    #[test]
    fn hv_roundtrip_file() {
        let idx = build_from_seqs(
            vec![(
                String::from("c"),
                b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            )],
            10,
            4,
            2,
            true,
            false,
        )
        .unwrap();
        let hv = index_to_hv(&idx, 1024, 3);
        let mut buf = Vec::new();
        write_hv(&mut buf, "test", 10, 1024, 3, idx.n_unique() as usize, &hv).unwrap();
        let loaded = read_hv(&mut std::io::Cursor::new(&buf)).unwrap();
        assert_eq!(loaded.name, "test");
        assert_eq!(loaded.k, 10);
        assert_eq!(loaded.dim, 1024);
        assert_eq!(loaded.sparse, 3);
        assert_eq!(loaded.n_kmer, idx.n_unique() as usize);
        assert_eq!(loaded.hv, hv);
    }

    #[test]
    fn crafted_hv_header_rejected_not_panic() {
        // Regression: dim = u32::MAX used to allocate ~17 GiB (abort) via
        // `vec![0i32; dim]`; crafted headers must error instead.
        let mut buf = Vec::new();
        buf.extend_from_slice(HV_MAGIC);
        buf.extend_from_slice(&HV_VERSION.to_le_bytes());
        buf.extend_from_slice(&10u32.to_le_bytes()); // k
        buf.extend_from_slice(&u32::MAX.to_le_bytes()); // dim
        buf.extend_from_slice(&3u32.to_le_bytes()); // sparse
        buf.extend_from_slice(&100u64.to_le_bytes()); // n_kmer
        buf.extend_from_slice(&4u32.to_le_bytes()); // name length
        buf.extend_from_slice(b"test");
        let err = read_hv(&mut std::io::Cursor::new(&buf))
            .err()
            .expect("crafted hv header must fail");
        // With memory overcommit the 17 GiB reserve may succeed and the
        // subsequent read fails on the truncated body; without it the reserve
        // itself errors. Either way it must be a friendly error, not an abort.
        let msg = err.to_string();
        assert!(
            msg.contains("hv dimension too large") || msg.contains("failed to fill whole buffer"),
            "got: {msg}"
        );
    }
}
