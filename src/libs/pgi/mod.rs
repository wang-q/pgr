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
    /// Per-key grouped, in entry order: packed `(contig_id, pos, strand)`
    /// records (see `pack_position`).
    pub positions: Vec<u64>,
}

/// Index header: parameters plus the contig table, without the k-mer records.
#[derive(Debug, Clone)]
pub struct PgiHeader {
    /// K-mer length (bp).
    pub k: usize,
    /// Syncmer length (bp).
    pub smer: usize,
    /// Syncmer window (bp).
    pub window: usize,
    /// `(name, length)` pairs in file order.
    pub contigs: Vec<(String, u64)>,
}

/// Per-record byte layout of the v2 packed occurrence stream.
#[derive(Debug, Clone, Copy)]
struct RecordLayout {
    kmer_bytes: usize,
    pos_bytes: usize,
    cont_bytes: usize,
}

impl RecordLayout {
    fn size(self) -> usize {
        self.kmer_bytes + self.pos_bytes + self.cont_bytes
    }
}

/// Read and validate the header (magic/version/params/contigs) plus the
/// per-record layout; returns the number of occurrence records.
fn read_header<R: Read>(r: &mut R) -> anyhow::Result<(PgiHeader, usize, RecordLayout)> {
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

    let mut contigs = Vec::with_capacity(n_contigs);
    for _ in 0..n_contigs {
        let nb = read_u32(r)? as usize;
        let mut name = vec![0u8; nb];
        r.read_exact(&mut name).context("reading contig name")?;
        let name = String::from_utf8(name).context("contig name utf8")?;
        let len = read_u64(r)?;
        contigs.push((name, len));
    }

    Ok((
        PgiHeader {
            k,
            smer,
            window,
            contigs,
        },
        n_records,
        RecordLayout {
            kmer_bytes,
            pos_bytes,
            cont_bytes,
        },
    ))
}

/// Decode one packed occurrence record into `(kmer, contig_id, pos, strand)`.
fn parse_record(rec: &[u8], k: usize, layout: RecordLayout) -> (u128, u32, u32, u8) {
    let kmer = unpack_kmer(&rec[..layout.kmer_bytes], k);
    let mut pos: u32 = 0;
    for (i, byte) in rec[layout.kmer_bytes..layout.kmer_bytes + layout.pos_bytes]
        .iter()
        .enumerate()
    {
        pos |= (*byte as u32) << (8 * i);
    }
    let mut cont: u32 = 0;
    for (i, byte) in rec[layout.kmer_bytes + layout.pos_bytes..]
        .iter()
        .enumerate()
    {
        cont |= (*byte as u32) << (8 * i);
    }
    let strand_bit = 1u32 << (8 * layout.cont_bytes - 1);
    let strand = (cont & strand_bit != 0) as u8;
    let cid = cont & (strand_bit - 1);
    (kmer, cid, pos, strand)
}

/// Bit width of the contig id within a packed position record (2^20 contigs).
const CID_BITS: u32 = 20;
/// Mask isolating the position bits of a packed position record.
const POS_MASK: u64 = (1 << 32) - 1;
/// Bit offset of the strand flag within a packed position record.
const STRAND_OFF: u32 = 32 + CID_BITS;

/// Pack `(contig_id, pos, strand)` into a `u64` position record.
fn pack_position(cid: u32, pos: u32, strand: u8) -> u64 {
    debug_assert!(cid >> CID_BITS == 0, "contig id exceeds packed range");
    (pos as u64) | ((cid as u64) << 32) | (((strand & 1) as u64) << STRAND_OFF)
}

/// Unpack a packed `u64` position record into `(contig_id, pos, strand)`.
fn unpack_position(rec: u64) -> (u32, u32, u8) {
    (
        ((rec >> 32) & ((1 << CID_BITS) - 1)) as u32,
        (rec & POS_MASK) as u32,
        ((rec >> STRAND_OFF) & 1) as u8,
    )
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
            for &prec in &self.positions[e.pos_start as usize..end] {
                let (cid, pos, strand) = unpack_position(prec);
                pack_kmer(e.kmer, self.k, &mut rec[..kmer_bytes]);
                let pb = &mut rec[kmer_bytes..kmer_bytes + pos_bytes];
                pb.fill(0);
                for (i, byte) in pb.iter_mut().enumerate() {
                    *byte = ((pos >> (8 * i)) & 0xff) as u8;
                }
                let cont = cid | ((strand as u32 & 1) << (8 * cont_bytes - 1));
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
        let (header, n_records, layout) = read_header(r)?;
        let PgiHeader {
            k,
            smer,
            window,
            contigs,
        } = header;
        let rec_size = layout.size();
        let mut entries: Vec<PgiEntry> = Vec::new();
        let mut positions: Vec<u64> = Vec::with_capacity(n_records);
        let mut last_kmer: Option<u128> = None;
        // Read the records in large chunks and parse from the slice (a
        // per-record `read_exact` through the trait object costs a virtual
        // dispatch for every one of the millions of records).
        let mut recs_left = n_records;
        let mut buf: Vec<u8> = Vec::with_capacity(1 << 20);
        while recs_left > 0 {
            let want = (recs_left * rec_size).min(buf.capacity()) / rec_size * rec_size;
            anyhow::ensure!(want > 0, "truncated index records");
            buf.resize(want, 0);
            r.read_exact(&mut buf).context("reading index records")?;
            for rec in buf.chunks_exact(rec_size) {
                let (kmer, cid, pos, strand) = parse_record(rec, k, layout);
                let pos_start = positions.len() as u32;
                positions.push(pack_position(cid, pos, strand));
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
                recs_left -= 1;
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

/// Streaming reader over a `.pgi` record stream: yields batches of complete
/// entries (each with its positions) without materializing the whole index.
/// The merge scans the reference index sequentially, so only the query index
/// needs to be fully resident. `pos_start` in yielded entries is 0 (the
/// positions ride along in the batch).
pub struct PgiStream<R: Read> {
    reader: R,
    header: PgiHeader,
    layout: RecordLayout,
    recs_left: usize,
    buf: Vec<u8>,
    buf_off: usize,
    buf_end: usize,
    pending: Option<(PgiEntry, Vec<u64>)>,
}

/// Stream read buffer size (1 MiB of records).
const STREAM_BUF: usize = 1 << 20;

impl<R: Read> PgiStream<R> {
    /// Open a stream over a .pgi reader, validating the header.
    pub fn open(mut reader: R) -> anyhow::Result<Self> {
        let (header, recs_left, layout) = read_header(&mut reader)?;
        Ok(Self {
            reader,
            header,
            layout,
            recs_left,
            buf: Vec::with_capacity(STREAM_BUF),
            buf_off: 0,
            buf_end: 0,
            pending: None,
        })
    }

    /// Index parameters and contig table.
    pub fn header(&self) -> &PgiHeader {
        &self.header
    }

    /// Next batch of up to `max_entries` complete entries (with positions),
    /// or an empty vec at end of stream. Entries are never split across
    /// batches.
    pub fn next_batch(&mut self, max_entries: usize) -> anyhow::Result<Vec<(PgiEntry, Vec<u64>)>> {
        let mut out: Vec<(PgiEntry, Vec<u64>)> = Vec::new();
        let mut cur: Option<(PgiEntry, Vec<u64>)> = self.pending.take();
        while self.recs_left > 0 && out.len() < max_entries {
            if self.buf_off == self.buf_end {
                let want = (self.recs_left * self.layout.size()).min(STREAM_BUF)
                    / self.layout.size()
                    * self.layout.size();
                anyhow::ensure!(want > 0, "truncated index records");
                self.buf.resize(want, 0);
                self.reader
                    .read_exact(&mut self.buf)
                    .context("reading index records")?;
                self.buf_off = 0;
                self.buf_end = want;
            }
            let rec = &self.buf[self.buf_off..self.buf_off + self.layout.size()];
            self.buf_off += self.layout.size();
            self.recs_left -= 1;
            let (kmer, cid, pos, strand) = parse_record(rec, self.header.k, self.layout);
            if let Some((e, poss)) = &mut cur {
                if e.kmer == kmer {
                    poss.push(pack_position(cid, pos, strand));
                    continue;
                }
                let (mut e, poss) = cur.take().expect("cur is Some");
                e.freq = poss.len() as u32;
                out.push((e, poss));
            }
            cur = Some((
                PgiEntry {
                    kmer,
                    pos_start: 0,
                    freq: 0,
                },
                vec![pack_position(cid, pos, strand)],
            ));
        }
        if self.recs_left == 0 {
            // End of stream: the carried entry is complete, yield it.
            if let Some((mut e, poss)) = cur {
                e.freq = poss.len() as u32;
                out.push((e, poss));
            }
        } else {
            self.pending = cur;
        }
        Ok(out)
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
    fn pack_unpack_position_roundtrip() {
        assert_eq!(unpack_position(pack_position(0, 0, 0)), (0, 0, 0));
        assert_eq!(
            unpack_position(pack_position(7, 4_641_652, 1)),
            (7, 4_641_652, 1)
        );
        let max_cid = (1 << CID_BITS) - 1;
        let rec = pack_position(max_cid, u32::MAX, 1);
        assert_eq!(unpack_position(rec), (max_cid, u32::MAX, 1));
        assert_eq!(pack_position(0, 0, 0), 0);
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

    #[test]
    fn stream_matches_full_read() {
        let seq: Vec<u8> = (0..200u32).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        let idx = build_from_seqs(vec![(String::from("c1"), seq)], 10, 4, 2, false).unwrap();
        let mut buf = Vec::new();
        idx.write(&mut buf).unwrap();
        let loaded = PgiIndex::read(&mut std::io::Cursor::new(&buf)).unwrap();

        let mut stream = PgiStream::open(std::io::Cursor::new(&buf)).unwrap();
        let mut entries = Vec::new();
        let mut positions = Vec::new();
        loop {
            let batch = stream.next_batch(64).unwrap();
            if batch.is_empty() {
                break;
            }
            for (e, poss) in batch {
                let start = positions.len() as u32;
                entries.push(PgiEntry {
                    pos_start: start,
                    ..e
                });
                positions.extend(poss);
            }
        }
        assert_eq!(stream.header().k, idx.k);
        assert_eq!(stream.header().contigs, idx.contigs);
        assert_eq!(entries, loaded.entries);
        assert_eq!(positions, loaded.positions);
    }
}
