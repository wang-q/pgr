//! Memory-mapped, zero-copy `.pgi` index view (FastGA GIX model).

use super::{
    pack_kmer, pack_position, parse_header_bytes, parse_record, unpack_kmer, PgiHeader, PgiQuery,
    Positions, RecordLayout,
};
use anyhow::Context;
use memmap2::{Mmap, MmapOptions};
use std::fs::File;
use std::path::Path;

/// Zero-copy, memory-mapped view of a `.pgi` index.
///
/// Records stay in the mapped pages: entries are located by binary search
/// over the packed (big-endian) k-mer bytes and positions are decoded on
/// demand, so the query index in `pgr align pgi` never materializes its
/// position table in memory.
pub struct PgiMmap {
    map: Mmap,
    k: usize,
    smer: usize,
    window: usize,
    contigs: Vec<(String, u64)>,
    layout: RecordLayout,
    records_off: usize,
    n_records: usize,
}

impl PgiMmap {
    /// Open and validate an index file, mapping it read-only.
    ///
    /// Uses `map_copy_read_only` (MAP_PRIVATE); the caller must not modify
    /// or truncate the file while the mapping is alive. The align command
    /// opens the query index read-only and never writes to it, which is the
    /// documented safe usage.
    pub fn open(path: &Path) -> anyhow::Result<Self> {
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        // SAFETY: the mapped file is opened read-only and is never modified
        // or truncated while this mapping exists (the caller only reads the
        // index); page faults cannot race with external writes to the file.
        let map = unsafe {
            MmapOptions::new()
                .map_copy_read_only(&file)
                .with_context(|| format!("mmap {}", path.display()))?
        };
        let (header, n_records, layout, records_off) = parse_header_bytes(&map)?;
        let PgiHeader {
            k,
            smer,
            window,
            contigs,
        } = header;
        let recs_len = n_records
            .checked_mul(layout.size())
            .context("index record region overflow")?;
        let rec_end = records_off
            .checked_add(recs_len)
            .context("index record region overflow")?;
        anyhow::ensure!(map.len() >= rec_end, "truncated index records");
        Ok(Self {
            map,
            k,
            smer,
            window,
            contigs,
            layout,
            records_off,
            n_records,
        })
    }

    /// K-mer length (bp).
    pub fn k(&self) -> usize {
        self.k
    }

    /// Syncmer length (bp).
    pub fn smer(&self) -> usize {
        self.smer
    }

    /// Syncmer window (bp).
    pub fn window(&self) -> usize {
        self.window
    }

    /// `(name, length)` pairs in file order.
    pub fn contigs(&self) -> &[(String, u64)] {
        &self.contigs
    }

    /// Total number of k-mer occurrence records.
    pub fn n_records(&self) -> usize {
        self.n_records
    }

    /// Bytes of one on-disk occurrence record.
    fn rec_bytes(&self, rec: usize) -> &[u8] {
        let off = self.records_off + rec * self.layout.size();
        &self.map[off..off + self.layout.size()]
    }

    /// Big-endian packed k-mer bytes of one record.
    fn rec_kmer_bytes(&self, rec: usize) -> &[u8] {
        &self.rec_bytes(rec)[..self.layout.kmer_bytes]
    }

    /// K-mer of one record.
    fn rec_kmer(&self, rec: usize) -> u128 {
        unpack_kmer(self.rec_kmer_bytes(rec), self.k)
    }

    /// First record index whose k-mer is `>= key` (records are grouped and
    /// ascending by packed k-mer bytes).
    fn lower_bound(&self, key: u128) -> usize {
        // Prefix ranges can be "one past the top" (lo + r == 2^(2k)); the
        // packed field cannot represent that sentinel, so cap it at the
        // record count like a partition_point over all keys.
        if 2 * self.k < 128 && key >= (1u128 << (2 * self.k)) {
            return self.n_records;
        }
        let mut kb = [0u8; 16];
        pack_kmer(key, self.k, &mut kb);
        let kb = &kb[..self.layout.kmer_bytes];
        let (mut lo, mut hi) = (0usize, self.n_records);
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if self.rec_kmer_bytes(mid) < kb {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        lo
    }

    /// One past the last record sharing entry `start`'s k-mer.
    fn group_end(&self, start: usize) -> usize {
        let kb = self.rec_kmer_bytes(start);
        let mut i = start + 1;
        while i < self.n_records && self.rec_kmer_bytes(i) == kb {
            i += 1;
        }
        i
    }
}

impl PgiQuery for PgiMmap {
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
        (self.lower_bound(lo), self.lower_bound(hi))
    }

    fn entry_next(&self, i: usize) -> usize {
        self.group_end(i)
    }

    fn entry_kmer(&self, i: usize) -> u128 {
        self.rec_kmer(i)
    }

    fn entry_freq(&self, i: usize) -> u32 {
        (self.group_end(i) - i) as u32
    }

    fn entry_positions(&self, i: usize) -> Positions<'_> {
        Positions::Mmap(MmapPosIter {
            map: &self.map,
            records_off: self.records_off,
            rec: i,
            end: self.group_end(i),
            k: self.k,
            layout: self.layout,
        })
    }
}

/// On-demand decoder over one entry's occurrence records in a mapped index.
pub struct MmapPosIter<'a> {
    map: &'a [u8],
    records_off: usize,
    rec: usize,
    end: usize,
    k: usize,
    layout: RecordLayout,
}

impl Iterator for MmapPosIter<'_> {
    type Item = u64;

    fn next(&mut self) -> Option<u64> {
        if self.rec == self.end {
            return None;
        }
        let off = self.records_off + self.rec * self.layout.size();
        let rec = &self.map[off..off + self.layout.size()];
        self.rec += 1;
        let (_, cid, pos, strand) = parse_record(rec, self.k, self.layout);
        Some(pack_position(cid, pos, strand))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let n = self.end - self.rec;
        (n, Some(n))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::pgi::build::build_from_seqs;
    use crate::libs::pgi::PgiIndex;
    use std::io::Cursor;

    fn seq(n: u32, seed: u32) -> Vec<u8> {
        let mut x = seed as u64;
        (0..n)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize & 3]
            })
            .collect()
    }

    #[test]
    fn mmap_matches_full_read() {
        let idx = build_from_seqs(
            vec![
                (String::from("c1"), seq(2000, 1)),
                (String::from("c2"), seq(1500, 2)),
            ],
            10,
            4,
            2,
            false,
            false,
        )
        .unwrap();
        let mut buf = Vec::new();
        idx.write(&mut buf).unwrap();
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("idx.pgi");
        std::fs::write(&path, &buf).unwrap();

        let loaded = PgiIndex::read(&mut Cursor::new(&buf)).unwrap();
        let mapped = PgiMmap::open(&path).unwrap();
        assert_eq!(mapped.k(), loaded.k);
        assert_eq!(mapped.smer(), loaded.smer);
        assert_eq!(mapped.window(), loaded.window);
        assert_eq!(mapped.contigs(), loaded.contigs);
        assert_eq!(mapped.n_records(), loaded.n_positions() as usize);

        let (_, end) = mapped.entry_range(0, u128::MAX);
        let mut rec = 0usize;
        for (ei, e) in loaded.entries.iter().enumerate() {
            assert_eq!(mapped.entry_kmer(rec), e.kmer, "entry {ei}");
            assert_eq!(mapped.entry_freq(rec), e.freq, "entry {ei}");
            let poss: Vec<u64> = mapped.entry_positions(rec).collect();
            let exp = &loaded.positions[e.pos_start as usize..(e.pos_start + e.freq) as usize];
            assert_eq!(poss, exp, "entry {ei}");
            rec = mapped.group_end(rec);
        }
        assert_eq!(rec, end);
    }

    #[test]
    fn mmap_entry_range_partitions_by_key() {
        let idx = build_from_seqs(
            vec![(String::from("c"), seq(4000, 7))],
            12,
            4,
            2,
            false,
            false,
        )
        .unwrap();
        let mut buf = Vec::new();
        idx.write(&mut buf).unwrap();
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("idx.pgi");
        std::fs::write(&path, &buf).unwrap();

        let mapped = PgiMmap::open(&path).unwrap();
        // A full-width range covers every entry; a half-width range splits
        // the key space exactly as the resident index's partition_point does.
        let (lo, hi) = mapped.entry_range(0, u128::MAX);
        assert_eq!(lo, 0);
        assert_eq!(hi, mapped.n_records());
        for &(qlo, qhi) in &[(0u128, 1u128 << 23), (1u128 << 23, 1u128 << 24)] {
            let (m0, m1) = mapped.entry_range(qlo, qhi);
            let (r0, r1) = (
                idx.entries.partition_point(|e| e.kmer < qlo),
                idx.entries.partition_point(|e| e.kmer < qhi),
            );
            // Group starts in the mmap are the same k-mers as the resident
            // entries; compare decoded keys, treating an end boundary as
            // "no key".
            let m_key = |i: usize| {
                if i < mapped.n_records() {
                    Some(mapped.entry_kmer(i))
                } else {
                    None
                }
            };
            let r_key = |i: usize| idx.entries.get(i).map(|e| e.kmer);
            assert_eq!(m_key(m0), r_key(r0));
            assert_eq!(m_key(m1), r_key(r1));
        }
    }

    #[test]
    fn mmap_truncated_records_rejected() {
        let idx = build_from_seqs(
            vec![(String::from("c"), seq(2000, 3))],
            10,
            4,
            2,
            false,
            false,
        )
        .unwrap();
        let mut buf = Vec::new();
        idx.write(&mut buf).unwrap();
        buf.truncate(buf.len() / 2);
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("idx.pgi");
        std::fs::write(&path, &buf).unwrap();
        let err = match PgiMmap::open(&path) {
            Ok(_) => panic!("truncated index must be rejected"),
            Err(e) => e,
        };
        assert!(err.to_string().contains("truncated"), "got: {err}");
    }
}
