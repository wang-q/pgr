//! pgr genome index (.pgi): syncmer-sparse sorted k-mer index.
//!
//! A `.pgi` file stores a genome's syncmer-sampled k-mers as a sorted,
//! duplicate-free table with per-key positions, supporting two-index merges
//! (distance / seed discovery) and hypervector projection. Design notes:
//! `notes/design/pbit.md`.

pub mod align;
pub mod build;
pub mod dist;
pub mod mmap;
pub mod to_hv;

pub use mmap::PgiMmap;

use anyhow::Context;
use mmap::MmapPosIter;
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
pub(crate) struct RecordLayout {
    kmer_bytes: usize,
    pos_bytes: usize,
    cont_bytes: usize,
}

impl RecordLayout {
    fn size(self) -> usize {
        self.kmer_bytes + self.pos_bytes + self.cont_bytes
    }
}

/// Fixed-size header fields (magic .. reserved), before the contig table.
const HEADER_FIXED: usize = 48;

/// Bounds-checked slice take that advances `off`.
fn take_bytes<'a>(buf: &'a [u8], off: &mut usize, n: usize) -> anyhow::Result<&'a [u8]> {
    let end = off.checked_add(n).context("header offset overflow")?;
    let s = buf.get(*off..end).context("truncated pgi header")?;
    *off = end;
    Ok(s)
}

/// Parse and validate the header (magic/version/params/contigs) plus the
/// per-record layout from a byte slice; returns the header, the number of
/// occurrence records, the layout, and the bytes consumed by the header.
pub(crate) fn parse_header_bytes(
    buf: &[u8],
) -> anyhow::Result<(PgiHeader, usize, RecordLayout, usize)> {
    let mut off = 0usize;
    let magic = take_bytes(buf, &mut off, 4)?;
    if magic != PGI_MAGIC {
        anyhow::bail!("not a pgr genome index (bad magic)");
    }
    let version = u32::from_le_bytes(take_bytes(buf, &mut off, 4)?.try_into().unwrap());
    if version != PGI_VERSION {
        anyhow::bail!("unsupported pgi version {version} (expected {PGI_VERSION})");
    }
    let k = u32::from_le_bytes(take_bytes(buf, &mut off, 4)?.try_into().unwrap()) as usize;
    let smer = u32::from_le_bytes(take_bytes(buf, &mut off, 4)?.try_into().unwrap()) as usize;
    let window = u32::from_le_bytes(take_bytes(buf, &mut off, 4)?.try_into().unwrap()) as usize;
    let n_contigs = u32::from_le_bytes(take_bytes(buf, &mut off, 4)?.try_into().unwrap()) as usize;
    let n_records = u64::from_le_bytes(take_bytes(buf, &mut off, 8)?.try_into().unwrap()) as usize;
    let kmer_bytes = u32::from_le_bytes(take_bytes(buf, &mut off, 4)?.try_into().unwrap()) as usize;
    let pos_bytes = u32::from_le_bytes(take_bytes(buf, &mut off, 4)?.try_into().unwrap()) as usize;
    let cont_bytes = u32::from_le_bytes(take_bytes(buf, &mut off, 4)?.try_into().unwrap()) as usize;
    let _reserved = u32::from_le_bytes(take_bytes(buf, &mut off, 4)?.try_into().unwrap());
    anyhow::ensure!(
        kmer_bytes == k.div_ceil(4),
        "bad kmer_bytes {kmer_bytes} for k={k}"
    );
    anyhow::ensure!((1..=4).contains(&pos_bytes), "bad pos_bytes {pos_bytes}");
    anyhow::ensure!((1..=4).contains(&cont_bytes), "bad cont_bytes {cont_bytes}");
    // Self-built indexes are capped at u16::MAX contigs and u32::MAX records
    // (see `build_from_seqs`); rejecting implausible counts turns crafted
    // headers into errors instead of allocation overflow/aborts.
    anyhow::ensure!(
        n_contigs <= u16::MAX as usize,
        "implausible contig count {n_contigs}"
    );
    anyhow::ensure!(
        n_records <= u32::MAX as usize,
        "implausible record count {n_records}"
    );

    let mut contigs = Vec::with_capacity(n_contigs);
    for _ in 0..n_contigs {
        let nb = u32::from_le_bytes(take_bytes(buf, &mut off, 4)?.try_into().unwrap()) as usize;
        let name = take_bytes(buf, &mut off, nb)?;
        let name = String::from_utf8(name.to_vec()).context("contig name utf8")?;
        let len = u64::from_le_bytes(take_bytes(buf, &mut off, 8)?.try_into().unwrap());
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
        off,
    ))
}

/// Read and validate the header (magic/version/params/contigs) plus the
/// per-record layout; returns the number of occurrence records.
fn read_header<R: Read>(r: &mut R) -> anyhow::Result<(PgiHeader, usize, RecordLayout)> {
    let mut buf = vec![0u8; HEADER_FIXED];
    r.read_exact(&mut buf).context("reading header")?;
    let n_contigs = u32::from_le_bytes(buf[20..24].try_into().unwrap()) as usize;
    anyhow::ensure!(
        n_contigs <= u16::MAX as usize,
        "implausible contig count {n_contigs}"
    );
    buf.reserve(n_contigs * 16);
    for _ in 0..n_contigs {
        let mut nb_bytes = [0u8; 4];
        r.read_exact(&mut nb_bytes)
            .context("reading contig name length")?;
        let nb = u32::from_le_bytes(nb_bytes) as usize;
        buf.extend_from_slice(&nb_bytes);
        let start = buf.len();
        buf.resize(start + nb, 0);
        r.read_exact(&mut buf[start..])
            .context("reading contig name")?;
        let mut len = [0u8; 8];
        r.read_exact(&mut len).context("reading contig length")?;
        buf.extend_from_slice(&len);
    }
    parse_header_bytes(&buf).map(|(h, n, l, _)| (h, n, l))
}

/// Decode one packed occurrence record into `(kmer, contig_id, pos, strand)`.
pub(crate) fn parse_record(rec: &[u8], k: usize, layout: RecordLayout) -> (u128, u32, u32, u8) {
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

/// Validate one decoded occurrence record against the index contig table.
///
/// The contig id must exist and the k-mer must start within the contig.
/// Crafted indexes must fail with a friendly error here instead of
/// panicking on the later contig-table lookups (Zero Panic).
pub(crate) fn validate_record(
    cid: u32,
    pos: u32,
    k: usize,
    contigs: &[(String, u64)],
) -> anyhow::Result<()> {
    let Some((_, len)) = contigs.get(cid as usize) else {
        anyhow::bail!(
            "index record contig id {cid} out of range ({} contigs)",
            contigs.len()
        );
    };
    anyhow::ensure!(
        pos as u64 + k as u64 <= *len,
        "index record position {pos} beyond contig {cid} length {len}"
    );
    Ok(())
}

/// Bit width of the contig id within a packed position record (2^20 contigs).
const CID_BITS: u32 = 20;
/// Mask isolating the position bits of a packed position record.
const POS_MASK: u64 = (1 << 32) - 1;
/// Bit offset of the strand flag within a packed position record.
const STRAND_OFF: u32 = 32 + CID_BITS;

/// Pack `(contig_id, pos, strand)` into a `u64` position record.
pub(crate) fn pack_position(cid: u32, pos: u32, strand: u8) -> u64 {
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

/// Read-only view over a pgi index's k-mer table, shared by the resident
/// [`PgiIndex`] and the memory-mapped [`PgiMmap`] query path.
///
/// Entries are addressable by an opaque `usize` index returned from
/// [`PgiQuery::entry_range`]; resident and mapped indexes differ in how
/// positions are stored (materialized `u64` records vs packed pages), so
/// only the decoded accessors are exposed here.
pub trait PgiQuery {
    /// K-mer length (bp).
    fn k(&self) -> usize;
    /// Syncmer length (bp).
    fn smer(&self) -> usize;
    /// Syncmer window (bp).
    fn window(&self) -> usize;
    /// `(name, length)` pairs in file order.
    fn contigs(&self) -> &[(String, u64)];
    /// Entry range whose k-mers lie in `[lo, hi)`.
    fn entry_range(&self, lo: u128, hi: u128) -> (usize, usize);
    /// Index one past the entry at `i` (entries may have a non-unit stride
    /// in the underlying storage).
    fn entry_next(&self, i: usize) -> usize;
    /// K-mer of the entry at index `i`.
    fn entry_kmer(&self, i: usize) -> u128;
    /// Position count of the entry at index `i`.
    fn entry_freq(&self, i: usize) -> u32;
    /// Packed position records of the entry at index `i`.
    fn entry_positions(&self, i: usize) -> Positions<'_>;
}

/// Position records of one entry: a resident slice or a decoder over the
/// mapped pages of a [`PgiMmap`].
pub enum Positions<'a> {
    /// Already-materialized packed records.
    Slice(std::slice::Iter<'a, u64>),
    /// Records decoded on demand from a memory-mapped index.
    Mmap(MmapPosIter<'a>),
}

impl Iterator for Positions<'_> {
    type Item = u64;

    fn next(&mut self) -> Option<u64> {
        match self {
            Positions::Slice(it) => it.next().copied(),
            Positions::Mmap(it) => it.next(),
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        match self {
            Positions::Slice(it) => it.size_hint(),
            Positions::Mmap(it) => it.size_hint(),
        }
    }
}

impl PgiQuery for PgiIndex {
    fn k(&self) -> usize {
        self.k
    }

    fn smer(&self) -> usize {
        self.smer
    }

    fn window(&self) -> usize {
        self.window
    }

    fn contigs(&self) -> &[(String, u64)] {
        &self.contigs
    }

    fn entry_range(&self, lo: u128, hi: u128) -> (usize, usize) {
        let i0 = self.entries.partition_point(|e| e.kmer < lo);
        let i1 = self.entries.partition_point(|e| e.kmer < hi);
        (i0, i1)
    }

    fn entry_next(&self, i: usize) -> usize {
        i + 1
    }

    fn entry_kmer(&self, i: usize) -> u128 {
        self.entries[i].kmer
    }

    fn entry_freq(&self, i: usize) -> u32 {
        self.entries[i].freq
    }

    fn entry_positions(&self, i: usize) -> Positions<'_> {
        let e = &self.entries[i];
        Positions::Slice(
            self.positions[e.pos_start as usize..(e.pos_start + e.freq) as usize].iter(),
        )
    }
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
        let mut positions: Vec<u64> = Vec::new();
        positions
            .try_reserve_exact(n_records)
            .context("allocating index positions")?;
        entries
            .try_reserve_exact(n_records)
            .context("allocating index entries")?;
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
                validate_record(cid, pos, k, &contigs)?;
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
            validate_record(cid, pos, self.header.k, &self.header.contigs)?;
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
pub(crate) fn pack_kmer(kmer: u128, k: usize, out: &mut [u8]) {
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
pub(crate) fn unpack_kmer(bytes: &[u8], k: usize) -> u128 {
    let mut x: u128 = 0;
    for &b in bytes {
        x = (x << 8) | b as u128;
    }
    x >> (8 * bytes.len() - 2 * k)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::pgi::build::build_from_seqs;

    fn crafted_header(n_contigs: u32, n_records: u64) -> Vec<u8> {
        let mut h = vec![0u8; HEADER_FIXED];
        h[0..4].copy_from_slice(PGI_MAGIC);
        h[4..8].copy_from_slice(&PGI_VERSION.to_le_bytes());
        h[8..12].copy_from_slice(&40u32.to_le_bytes());
        h[12..16].copy_from_slice(&8u32.to_le_bytes());
        h[16..20].copy_from_slice(&5u32.to_le_bytes());
        h[20..24].copy_from_slice(&n_contigs.to_le_bytes());
        h[24..32].copy_from_slice(&n_records.to_le_bytes());
        h[32..36].copy_from_slice(&10u32.to_le_bytes()); // kmer_bytes for k=40
        h[36..40].copy_from_slice(&3u32.to_le_bytes()); // pos_bytes
        h[40..44].copy_from_slice(&1u32.to_le_bytes()); // cont_bytes
        h
    }

    #[test]
    fn crafted_record_count_rejected_not_panic() {
        // Regression: a header claiming u64::MAX records used to hit
        // `Vec::with_capacity` capacity overflow (panic); it must error.
        let err =
            PgiIndex::read(&mut std::io::Cursor::new(crafted_header(0, u64::MAX))).unwrap_err();
        assert!(
            err.to_string().contains("implausible record count"),
            "got: {err}"
        );
    }

    #[test]
    fn crafted_contig_count_rejected_not_panic() {
        // Regression: a header claiming u32::MAX contigs used to reserve
        // ~64 GiB in `read_header` (allocation abort on most machines).
        let err =
            PgiIndex::read(&mut std::io::Cursor::new(crafted_header(u32::MAX, 0))).unwrap_err();
        assert!(
            err.to_string().contains("implausible contig count"),
            "got: {err}"
        );
    }

    #[test]
    fn crafted_record_contig_rejected_not_panic() {
        // Regression: an occurrence record with a contig id beyond the
        // contig table used to panic in the alignment contig lookups; the
        // resident reader and the streaming reader must reject it with a
        // friendly error.
        let seq: Vec<u8> = (0..200u32).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        let idx = build_from_seqs(vec![(String::from("c1"), seq)], 10, 4, 2, false).unwrap();
        let mut buf = Vec::new();
        idx.write(&mut buf).unwrap();
        corrupt_first_record_contig(&mut buf);

        let err = PgiIndex::read(&mut std::io::Cursor::new(&buf)).unwrap_err();
        assert!(err.to_string().contains("out of range"), "got: {err}");

        let mut stream = PgiStream::open(std::io::Cursor::new(&buf)).unwrap();
        let err = stream.next_batch(16).unwrap_err();
        assert!(err.to_string().contains("out of range"), "got: {err}");
    }

    /// Set the contig id of the first occurrence record to 0x7f (127), which
    /// is beyond any single-contig index's table.
    fn corrupt_first_record_contig(buf: &mut [u8]) {
        let (_h, _n, layout, records_off) =
            parse_header_bytes(buf).expect("valid test index header");
        let cont_off = records_off + layout.kmer_bytes + layout.pos_bytes;
        buf[cont_off] = 0x7f;
    }

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
