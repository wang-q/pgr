//! pgr genome index (.pgi): syncmer-sparse sorted k-mer index.
//!
//! A `.pgi` file stores a genome's syncmer-sampled k-mers as a sorted,
//! duplicate-free table with per-key positions, supporting two-index merges
//! (distance / seed discovery) and hypervector projection. Design notes:
//! `notes/design/pbit.md`.

pub mod align;
pub mod build;
pub mod dist;
pub mod to_hv;

use anyhow::Context;
use std::io::{Read, Write};

/// File magic.
pub const PGI_MAGIC: &[u8; 4] = b"PGI1";
/// Format version (v2: GIX-style packed per-occurrence records; see
/// notes/benchmarks/bench-pgi-vs-gix-storage.md).
pub const PGI_VERSION: u32 = 2;

/// One unique k-mer entry; its positions live in `PgiIndex::positions`
/// starting at `pos_start` and spanning `freq` records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PgiEntry {
    /// 2-bit encoded k-mer (high bits first), `2*k` significant bits.
    pub kmer: u128,
    /// Start index into `PgiIndex::positions`.
    pub pos_start: u32,
    /// Number of positions for this k-mer.
    pub freq: u32,
}

/// In-memory .pgi index.
#[derive(Debug, Clone, Default)]
pub struct PgiIndex {
    pub k: usize,
    pub smer: usize,
    pub window: usize,
    pub contigs: Vec<(String, u64)>,
    /// Sorted ascending by `kmer`.
    pub entries: Vec<PgiEntry>,
    /// Per-key grouped, in entry order: `(contig_id, pos, strand)`.
    pub positions: Vec<(u32, u32, u8)>,
}

impl PgiIndex {
    /// Number of unique k-mers.
    pub fn n_unique(&self) -> u64 {
        self.entries.len() as u64
    }

    /// Total number of k-mer positions.
    pub fn n_positions(&self) -> u64 {
        self.positions.len() as u64
    }

    /// Serialize to a writer (binary; v2 packed per-occurrence records).
    ///
    /// On disk each record is `kmer_bytes + pos_bytes + cont_bytes`:
    /// the 2-bit k-mer packed big-endian, the position in minimal little-
    /// endian bytes, and `contig_id | (strand << (8*cont_bytes-1))`.
    /// Entries are re-grouped on read, so the in-memory structure stays
    /// unchanged.
    pub fn write<W: Write>(&self, w: &mut W) -> anyhow::Result<()> {
        anyhow::ensure!(!self.contigs.is_empty(), "index has no contigs");
        let kmer_bytes = self.k.div_ceil(4);
        let max_len = self.contigs.iter().map(|(_, l)| *l).max().unwrap_or(0);
        let mut pos_bytes = 1usize;
        while (1u64 << (8 * pos_bytes)) < max_len {
            pos_bytes += 1;
        }
        let mut cont_bytes = 1usize;
        while (1u64 << (8 * cont_bytes)) < 2 * self.contigs.len() as u64 {
            cont_bytes += 1;
        }
        w.write_all(PGI_MAGIC)?;
        w.write_all(&PGI_VERSION.to_le_bytes())?;
        w.write_all(&(self.k as u32).to_le_bytes())?;
        w.write_all(&(self.smer as u32).to_le_bytes())?;
        w.write_all(&(self.window as u32).to_le_bytes())?;
        w.write_all(&(self.contigs.len() as u32).to_le_bytes())?;
        w.write_all(&self.n_positions().to_le_bytes())?;
        w.write_all(&(kmer_bytes as u32).to_le_bytes())?;
        w.write_all(&(pos_bytes as u32).to_le_bytes())?;
        w.write_all(&(cont_bytes as u32).to_le_bytes())?;
        w.write_all(&0u32.to_le_bytes())?; // reserved

        for (name, len) in &self.contigs {
            let nb = name.len() as u32;
            w.write_all(&nb.to_le_bytes())?;
            w.write_all(name.as_bytes())?;
            w.write_all(&len.to_le_bytes())?;
        }
        let mut rec = vec![0u8; kmer_bytes + pos_bytes + cont_bytes];
        for e in &self.entries {
            let end = (e.pos_start + e.freq) as usize;
            for (cid, pos, strand) in &self.positions[e.pos_start as usize..end] {
                pack_kmer(e.kmer, self.k, &mut rec[..kmer_bytes]);
                let pb = &mut rec[kmer_bytes..kmer_bytes + pos_bytes];
                pb.fill(0);
                for (i, byte) in pb.iter_mut().enumerate() {
                    *byte = ((*pos >> (8 * i)) & 0xff) as u8;
                }
                let cont = *cid | ((*strand as u32 & 1) << (8 * cont_bytes - 1));
                let cb = &mut rec[kmer_bytes + pos_bytes..];
                cb.fill(0);
                for (i, byte) in cb.iter_mut().enumerate() {
                    *byte = ((cont >> (8 * i)) & 0xff) as u8;
                }
                w.write_all(&rec)?;
            }
        }
        Ok(())
    }

    /// Deserialize from a reader; validates magic/version and re-groups the
    /// packed per-occurrence records into entries + positions.
    pub fn read<R: Read>(r: &mut R) -> anyhow::Result<Self> {
        let mut magic = [0u8; 4];
        r.read_exact(&mut magic).context("reading magic")?;
        if &magic != PGI_MAGIC {
            anyhow::bail!("not a pgr genome index (bad magic)");
        }
        let version = read_u32(r)?;
        if version != PGI_VERSION {
            anyhow::bail!("unsupported pgi version {version} (expected {PGI_VERSION})");
        }
        let k = read_u32(r)? as usize;
        let smer = read_u32(r)? as usize;
        let window = read_u32(r)? as usize;
        let n_contigs = read_u32(r)? as usize;
        let n_records = read_u64(r)? as usize;
        let kmer_bytes = read_u32(r)? as usize;
        let pos_bytes = read_u32(r)? as usize;
        let cont_bytes = read_u32(r)? as usize;
        let _reserved = read_u32(r)?;
        anyhow::ensure!(
            kmer_bytes == k.div_ceil(4),
            "bad kmer_bytes {kmer_bytes} for k={k}"
        );
        anyhow::ensure!((1..=4).contains(&pos_bytes), "bad pos_bytes {pos_bytes}");
        anyhow::ensure!((1..=4).contains(&cont_bytes), "bad cont_bytes {cont_bytes}");
        let strand_bit = 1u32 << (8 * cont_bytes - 1);

        let mut contigs = Vec::with_capacity(n_contigs);
        for _ in 0..n_contigs {
            let nb = read_u32(r)? as usize;
            let mut name = vec![0u8; nb];
            r.read_exact(&mut name).context("reading contig name")?;
            let name = String::from_utf8(name).context("contig name utf8")?;
            let len = read_u64(r)?;
            contigs.push((name, len));
        }

        let mut rec = vec![0u8; kmer_bytes + pos_bytes + cont_bytes];
        let mut entries: Vec<PgiEntry> = Vec::new();
        let mut positions: Vec<(u32, u32, u8)> = Vec::with_capacity(n_records);
        let mut last_kmer: Option<u128> = None;
        for _ in 0..n_records {
            r.read_exact(&mut rec)?;
            let kmer = unpack_kmer(&rec[..kmer_bytes], k);
            let mut pos: u32 = 0;
            for (i, byte) in rec[kmer_bytes..kmer_bytes + pos_bytes].iter().enumerate() {
                pos |= (*byte as u32) << (8 * i);
            }
            let mut cont: u32 = 0;
            for (i, byte) in rec[kmer_bytes + pos_bytes..].iter().enumerate() {
                cont |= (*byte as u32) << (8 * i);
            }
            let strand = (cont & strand_bit != 0) as u8;
            let cid = cont & (strand_bit - 1);
            let pos_start = positions.len() as u32;
            positions.push((cid, pos, strand));
            if last_kmer == Some(kmer) {
                entries.last_mut().unwrap().freq += 1;
            } else {
                entries.push(PgiEntry {
                    kmer,
                    pos_start,
                    freq: 1,
                });
                last_kmer = Some(kmer);
            }
        }

        Ok(PgiIndex {
            k,
            smer,
            window,
            contigs,
            entries,
            positions,
        })
    }
}

/// Pack a 2-bit k-mer (high bits first) into `kmer_bytes = ceil(k/4)` bytes,
/// big-endian, low bits zero-padded.
fn pack_kmer(kmer: u128, k: usize, out: &mut [u8]) {
    for (i, byte) in out.iter_mut().enumerate() {
        // Bases 4i..4i+3 high-aligned within the byte; missing trailing
        // bases (k % 4 != 0) leave the low bits zero.
        let mut b = 0u8;
        for j in 0..4 {
            let base_idx = 4 * i + j;
            if base_idx < k {
                let base = ((kmer >> (2 * (k - 1 - base_idx))) & 3) as u8;
                b |= base << (2 * (3 - j));
            }
        }
        *byte = b;
    }
}

/// Unpack a big-endian 2-bit k-mer byte array back to a `u128`.
fn unpack_kmer(bytes: &[u8], k: usize) -> u128 {
    let mut x: u128 = 0;
    for &b in bytes {
        x = (x << 8) | b as u128;
    }
    x >> (8 * bytes.len() - 2 * k)
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
    fn pack_unpack_kmer_roundtrip() {
        for &k in &[10usize, 32, 40, 63, 64] {
            let nbytes = k.div_ceil(4);
            let mut x: u128 = 0x0123_4567_89ab_cdef_0123_4567_89ab_cdef;
            if 2 * k < 128 {
                x &= (1u128 << (2 * k)) - 1;
            }
            let mut buf = vec![0u8; nbytes];
            pack_kmer(x, k, &mut buf);
            assert_eq!(unpack_kmer(&buf, k), x, "k={k}");
        }
    }

    #[test]
    fn write_read_roundtrip() {
        let seq: Vec<u8> = (0..200u32).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        let idx = build_from_seqs(vec![(String::from("c1"), seq)], 10, 4, 2, false).unwrap();
        let mut buf = Vec::new();
        idx.write(&mut buf).unwrap();
        let loaded = PgiIndex::read(&mut std::io::Cursor::new(&buf)).unwrap();
        assert_eq!(loaded.k, idx.k);
        assert_eq!(loaded.contigs, idx.contigs);
        assert_eq!(loaded.entries, idx.entries);
        assert_eq!(loaded.positions, idx.positions);
    }
}
