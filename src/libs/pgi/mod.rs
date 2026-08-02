//! pgr genome index (.pgi): syncmer-sparse sorted k-mer index.
//!
//! A `.pgi` file stores a genome's syncmer-sampled k-mers as a sorted,
//! duplicate-free table with per-key positions, supporting two-index merges
//! (distance / seed discovery) and hypervector projection. Design notes:
//! `notes/design/pbit-index-extension.md`.

pub mod build;
pub mod dist;
pub mod to_hv;

use anyhow::Context;
use std::io::{Read, Write};

/// File magic.
pub const PGI_MAGIC: &[u8; 4] = b"PGI1";
/// Format version.
pub const PGI_VERSION: u32 = 1;

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

    /// Serialize to a writer (binary, little-endian).
    pub fn write<W: Write>(&self, w: &mut W) -> anyhow::Result<()> {
        w.write_all(PGI_MAGIC)?;
        w.write_all(&PGI_VERSION.to_le_bytes())?;
        w.write_all(&(self.k as u32).to_le_bytes())?;
        w.write_all(&(self.smer as u32).to_le_bytes())?;
        w.write_all(&(self.window as u32).to_le_bytes())?;
        w.write_all(&(self.contigs.len() as u32).to_le_bytes())?;
        w.write_all(&self.n_unique().to_le_bytes())?;
        w.write_all(&self.n_positions().to_le_bytes())?;
        w.write_all(&0u64.to_le_bytes())?; // reserved

        for (name, len) in &self.contigs {
            let nb = name.len() as u32;
            w.write_all(&nb.to_le_bytes())?;
            w.write_all(name.as_bytes())?;
            w.write_all(&len.to_le_bytes())?;
        }
        // Batch each fixed-size record into one write to cut call overhead.
        let mut rec = [0u8; 24];
        for e in &self.entries {
            rec[..16].copy_from_slice(&e.kmer.to_le_bytes());
            rec[16..20].copy_from_slice(&e.pos_start.to_le_bytes());
            rec[20..24].copy_from_slice(&e.freq.to_le_bytes());
            w.write_all(&rec)?;
        }
        let mut posrec = [0u8; 9];
        for (cid, pos, strand) in &self.positions {
            posrec[..4].copy_from_slice(&cid.to_le_bytes());
            posrec[4..8].copy_from_slice(&pos.to_le_bytes());
            posrec[8] = *strand;
            w.write_all(&posrec)?;
        }
        Ok(())
    }

    /// Deserialize from a reader; validates magic/version.
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
        let n_unique = read_u64(r)? as usize;
        let n_positions = read_u64(r)? as usize;
        let _reserved = read_u64(r)?;

        let mut contigs = Vec::with_capacity(n_contigs);
        for _ in 0..n_contigs {
            let nb = read_u32(r)? as usize;
            let mut name = vec![0u8; nb];
            r.read_exact(&mut name).context("reading contig name")?;
            let name = String::from_utf8(name).context("contig name utf8")?;
            let len = read_u64(r)?;
            contigs.push((name, len));
        }

        let mut entries = Vec::with_capacity(n_unique);
        for _ in 0..n_unique {
            let mut kb = [0u8; 16];
            r.read_exact(&mut kb)?;
            let kmer = u128::from_le_bytes(kb);
            let pos_start = read_u32(r)?;
            let freq = read_u32(r)?;
            entries.push(PgiEntry {
                kmer,
                pos_start,
                freq,
            });
        }

        let mut positions = Vec::with_capacity(n_positions);
        for _ in 0..n_positions {
            let cid = read_u32(r)?;
            let pos = read_u32(r)?;
            let mut strand = [0u8; 1];
            r.read_exact(&mut strand)?;
            positions.push((cid, pos, strand[0]));
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
