//! Decompressor: reads a `.pbit` archive and extracts reference / sample
//! sequences.
//!
//! Holds an `R: Read + Seek` reader directly (no archive wrapper). Parses the
//! header + footer + reference index + delta data headers + sample index at
//! construction time, then serves random-access queries via `get_contig` /
//! `get_sample` / `SequenceReader`.

use anyhow::{anyhow, Context, Result};
use indexmap::IndexMap;
use lru::LruCache;
use std::collections::HashSet;
use std::io::{Read, Seek, SeekFrom, Write};
use std::num::NonZeroUsize;
use std::path::Path;

use crate::libs::fmt::twobit::read_2bit_record;
use crate::libs::io::SequenceReader;
use crate::libs::nt;

use super::collection::Collection;
use super::format::{
    read_ref_index, read_u32_le, DeltaMeta, PbitFooter, PbitHeader, RefGroupEntry,
};
use super::segment::Segment;

/// Decompressor for a `.pbit` archive.
pub struct Decompressor<R: Read + Seek> {
    reader: R,
    header: PbitHeader,
    #[allow(dead_code)]
    footer: PbitFooter,
    ref_groups: Vec<RefGroupEntry>,
    /// contig name → Vec<ref_group_id> (reference segments, ordered).
    contig_groups: IndexMap<String, Vec<u32>>,
    /// All contig names appearing in any sample's collection (for
    /// `contains_contig`).
    contig_set: HashSet<String>,
    collection: Collection,
    /// delta_meta[ref_group_id][delta_id] → header info (no packed data).
    /// Kept for future `stat --deltas` exposure; decoding re-reads headers
    /// from disk via `delta_offsets`.
    #[allow(dead_code)]
    delta_meta: Vec<Vec<DeltaMeta>>,
    /// delta_offsets[ref_group_id][delta_id] → file offset of the delta's
    /// 9-byte header (followed by `packed_size` bytes).
    delta_offsets: Vec<Vec<u64>>,
    /// LRU cache: ref_group_id → decoded reference segment DNA (ASCII).
    ref_cache: LruCache<u32, Vec<u8>>,
    /// LRU cache: (ref_group_id, delta_id) → decoded raw sample segment.
    delta_cache: LruCache<(u32, u32), Vec<u8>>,
    min_match_len: u32,
}

impl Decompressor<std::io::BufReader<std::fs::File>> {
    /// Open from a file path (mirrors `TwoBitFile::open`).
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = std::fs::File::open(&path)
            .with_context(|| format!("failed to open pbit file: {}", path.as_ref().display()))?;
        let reader = std::io::BufReader::new(file);
        Self::new(reader)
    }
}

impl<R: Read + Seek> Decompressor<R> {
    /// Construct from an already-opened reader: parse header + footer +
    /// indexes + scan delta data headers.
    pub fn new(mut reader: R) -> Result<Self> {
        // Read header.
        let header = PbitHeader::read_from(&mut reader)?;

        // Read footer.
        let footer = PbitFooter::read_at_end(&mut reader)?;

        // Read reference index.
        reader.seek(SeekFrom::Start(footer.ref_index_offset))?;
        let ref_groups = read_ref_index(&mut reader)?;
        if ref_groups.len() as u32 != header.ref_group_count {
            return Err(anyhow!(
                "ref_group_count mismatch: header={}, index={}",
                header.ref_group_count,
                ref_groups.len()
            ));
        }

        // Build contig_groups (contig name → ref_group_ids in order).
        let mut contig_groups: IndexMap<String, Vec<u32>> = IndexMap::new();
        for (i, entry) in ref_groups.iter().enumerate() {
            contig_groups
                .entry(entry.contig_name.clone())
                .or_default()
                .push(i as u32);
        }

        // Scan delta data: read each delta's 9-byte header, build delta_meta
        // and delta_offsets (without decompressing data).
        reader.seek(SeekFrom::Start(footer.delta_data_offset))?;
        let ref_group_count = read_u32_le(&mut reader)? as usize;
        if ref_group_count != header.ref_group_count as usize {
            return Err(anyhow!(
                "ref_group_count mismatch: header={}, delta_data={}",
                header.ref_group_count,
                ref_group_count
            ));
        }
        let mut delta_meta: Vec<Vec<DeltaMeta>> = Vec::with_capacity(ref_group_count);
        let mut delta_offsets: Vec<Vec<u64>> = Vec::with_capacity(ref_group_count);
        for _ in 0..ref_group_count {
            let delta_count = read_u32_le(&mut reader)? as usize;
            let mut metas = Vec::with_capacity(delta_count);
            let mut offsets = Vec::with_capacity(delta_count);
            for _ in 0..delta_count {
                let offset = reader.stream_position()?;
                let meta = DeltaMeta::read_header(&mut reader)?;
                metas.push(meta);
                offsets.push(offset);
                // Skip the packed data.
                reader.seek(SeekFrom::Current(meta.packed_size as i64))?;
            }
            delta_meta.push(metas);
            delta_offsets.push(offsets);
        }

        // Read sample index (collection, flate2-compressed).
        reader.seek(SeekFrom::Start(footer.sample_index_offset))?;
        let mut compressed = Vec::new();
        reader.read_to_end(&mut compressed)?;
        // The footer is 24 bytes at the end; the collection is everything from
        // sample_index_offset to the footer. But we already read_to_end which
        // includes the footer. We need to strip the trailing 24 bytes.
        // Actually, the footer is already parsed via read_at_end, and the
        // collection bytes are between sample_index_offset and footer position.
        // Since we used read_to_end, we have collection + footer. Strip 24.
        let footer_len = 24;
        if compressed.len() > footer_len {
            compressed.truncate(compressed.len() - footer_len);
        }
        let collection = Collection::deserialize(&compressed)?;

        // Build contig_set from collection.
        let mut contig_set: HashSet<String> = HashSet::new();
        for (_sample, contigs) in &collection.samples {
            for cs in contigs {
                contig_set.insert(cs.contig_name.clone());
            }
        }

        let min_match_len = if header.kmer_len >= 4 {
            header.kmer_len + 3
        } else {
            18
        };

        Ok(Self {
            reader,
            header,
            footer,
            ref_groups,
            contig_groups,
            contig_set,
            collection,
            delta_meta,
            delta_offsets,
            ref_cache: LruCache::new(NonZeroUsize::new(64).unwrap()),
            delta_cache: LruCache::new(NonZeroUsize::new(256).unwrap()),
            min_match_len,
        })
    }

    /// Check if a contig name exists in any sample's collection.
    pub fn contains_contig(&self, name: &str) -> bool {
        self.contig_set.contains(name)
    }

    /// List all sample names.
    pub fn list_samples(&self) -> Vec<&str> {
        self.collection.list_samples()
    }

    /// List contig names for a sample (or all contigs across all samples if
    /// `sample` is `None`).
    pub fn list_contigs(&self, sample: Option<&str>) -> Vec<&str> {
        match sample {
            Some(s) => self.collection.list_contigs(s),
            None => {
                let mut seen: HashSet<&str> = HashSet::new();
                let mut out = Vec::new();
                for contigs in self.collection.samples.values() {
                    for cs in contigs {
                        if seen.insert(cs.contig_name.as_str()) {
                            out.push(cs.contig_name.as_str());
                        }
                    }
                }
                out
            }
        }
    }

    /// Return the reference group entries (for `stat --refs`).
    pub fn ref_groups(&self) -> &[RefGroupEntry] {
        &self.ref_groups
    }

    /// Return the header (for `stat` overview).
    pub fn header(&self) -> &PbitHeader {
        &self.header
    }

    /// Return the collection (for `stat --contigs`).
    pub fn collection(&self) -> &Collection {
        &self.collection
    }

    /// Read and decode a reference segment (2bit record) by ref_group_id.
    /// Uses an LRU cache to avoid re-reading.
    fn read_ref_segment(&mut self, ref_group_id: u32) -> Result<Vec<u8>> {
        if let Some(cached) = self.ref_cache.get(&ref_group_id) {
            return Ok(cached.clone());
        }
        let offset = self.ref_groups[ref_group_id as usize].segment_offset;
        self.reader.seek(SeekFrom::Start(offset))?;
        // Read the full reference segment (no slice, no masking removal).
        let seq = read_2bit_record(&mut self.reader, false, None, None, true)?;
        let seq_bytes = seq.into_bytes();
        self.ref_cache.put(ref_group_id, seq_bytes.clone());
        Ok(seq_bytes)
    }

    /// Read and flate2-decompress a delta's packed data, then LZ-diff decode
    /// it against the reference segment. Uses an LRU cache.
    fn decode_delta(&mut self, ref_group_id: u32, delta_id: u32) -> Result<Vec<u8>> {
        let key = (ref_group_id, delta_id);
        if let Some(cached) = self.delta_cache.get(&key) {
            return Ok(cached.clone());
        }

        // Read reference segment.
        let ref_dna = self.read_ref_segment(ref_group_id)?;

        // Read and decompress delta.
        let offset = self.delta_offsets[ref_group_id as usize][delta_id as usize];
        self.reader.seek(SeekFrom::Start(offset))?;
        let meta = DeltaMeta::read_header(&mut self.reader)?;
        let mut packed = vec![0u8; meta.packed_size as usize];
        self.reader.read_exact(&mut packed)?;
        let mut decoder = flate2::read::GzDecoder::new(&packed[..]);
        let mut delta = Vec::new();
        decoder.read_to_end(&mut delta)?;

        // LZ-diff decode.
        let mut seg = Segment::new(self.min_match_len);
        seg.prepare(&ref_dna);
        let mut decoded = seg.get(&delta)?;

        // Apply reverse-complement if needed.
        if meta.is_rev_comp {
            decoded = nt::rev_comp(&decoded).collect();
        }

        self.delta_cache.put(key, decoded.clone());
        Ok(decoded)
    }

    /// Extract a contig from ALL samples (getctg semantics), optionally sliced
    /// to `[start, end)`. Writes one FASTA entry per sample that has this
    /// contig. If `strand` is `"-"`, each sequence is reverse-complemented
    /// before writing.
    pub fn get_contig(
        &mut self,
        contig: &str,
        start: Option<usize>,
        end: Option<usize>,
        strand: &str,
        out: &mut impl Write,
    ) -> Result<()> {
        let line_width = 60;

        // Collect (sample_name, segments) pairs first to release the immutable
        // borrow on self.collection before calling self.decode_delta(). We
        // must clone the segments (not just borrow) to fully release the
        // immutable borrow on self.collection.
        let sample_segs: Vec<(String, Vec<crate::libs::pbit::collection::SegmentDesc>)> = self
            .collection
            .samples
            .iter()
            .filter_map(|(s, contigs)| {
                contigs
                    .iter()
                    .find(|c| c.contig_name == contig)
                    .map(|cs| (s.clone(), cs.segments.clone()))
            })
            .collect();

        for (sample, segments) in sample_segs {
            // Decode and concatenate all segments of this contig.
            let mut full_seq = Vec::new();
            for seg in &segments {
                let decoded = self.decode_delta(seg.ref_group_id, seg.delta_id)?;
                full_seq.extend_from_slice(&decoded);
            }

            // Apply slice [start, end).
            let total_len = full_seq.len();
            let s = start.unwrap_or(0).min(total_len);
            let e = end.unwrap_or(total_len).min(total_len);
            if s < e {
                let slice = &full_seq[s..e];
                let seq_bytes: Vec<u8> = if strand == "-" {
                    nt::rev_comp(slice).collect()
                } else {
                    slice.to_vec()
                };

                // Write FASTA header.
                let header = match (start, end, strand) {
                    (Some(_), Some(_), _) => {
                        format!(">{} {}:{}-{}({})", sample, contig, s + 1, e, strand)
                    }
                    _ => format!(">{} {}", sample, contig),
                };
                writeln!(out, "{}", header)?;
                write_fasta_seq(out, &seq_bytes, line_width)?;
            }
        }
        Ok(())
    }

    /// Extract all contigs of a single sample, writing FASTA entries.
    pub fn get_sample(&mut self, sample: &str, out: &mut impl Write) -> Result<()> {
        let line_width = 60;

        // Collect (contig_name, segments) pairs first to release the immutable
        // borrow on self.collection before calling self.decode_delta().
        let contig_segs: Vec<(String, Vec<crate::libs::pbit::collection::SegmentDesc>)> =
            match self.collection.samples.get(sample) {
                Some(c) => c
                    .iter()
                    .map(|cs| (cs.contig_name.clone(), cs.segments.clone()))
                    .collect(),
                None => {
                    return Err(anyhow!("sample '{}' not found in archive", sample));
                }
            };

        for (contig_name, segments) in contig_segs {
            let mut full_seq = Vec::new();
            for seg in &segments {
                let decoded = self.decode_delta(seg.ref_group_id, seg.delta_id)?;
                full_seq.extend_from_slice(&decoded);
            }
            writeln!(out, ">{}", contig_name)?;
            write_fasta_seq(out, &full_seq, line_width)?;
        }
        Ok(())
    }
}

/// Write a sequence byte slice as FASTA with line wrapping.
fn write_fasta_seq(out: &mut impl Write, seq: &[u8], line_width: usize) -> Result<()> {
    if line_width == 0 {
        out.write_all(seq)?;
        writeln!(out)?;
    } else {
        for chunk in seq.chunks(line_width) {
            out.write_all(chunk)?;
            writeln!(out)?;
        }
    }
    Ok(())
}

impl<R: Read + Seek> SequenceReader for Decompressor<R> {
    /// Read `[start, end)` from the reference sequence `name`. `None` means
    /// "from start" / "to end". Reads the REFERENCE layer (not sample layer).
    fn read_sequence(
        &mut self,
        name: &str,
        start: Option<usize>,
        end: Option<usize>,
    ) -> Result<String> {
        // Clone the ref_group_ids first to release the immutable borrow on
        // self.contig_groups before calling self.read_ref_segment().
        let ref_group_ids: Vec<u32> = match self.contig_groups.get(name) {
            Some(ids) => ids.clone(),
            None => return Err(anyhow!("contig '{}' not found in reference", name)),
        };

        // Walk segments, accumulate lengths, read only segments overlapping
        // [start, end).
        let mut result = Vec::new();
        let mut offset: usize = 0;
        let s = start.unwrap_or(0);
        let e = end.unwrap_or(usize::MAX);

        for rgid in ref_group_ids {
            let seg_dna = self.read_ref_segment(rgid)?;
            let seg_len = seg_dna.len();
            let seg_end = offset + seg_len;
            if seg_end > s && offset < e {
                // This segment overlaps [s, e).
                let local_start = s.saturating_sub(offset);
                let local_end = (e - offset).min(seg_len);
                result.extend_from_slice(&seg_dna[local_start..local_end]);
            }
            offset = seg_end;
            if offset >= e {
                break;
            }
        }

        Ok(String::from_utf8(result)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::pbit::compressor::Compressor;

    fn random_dna(len: usize, seed: u64) -> Vec<u8> {
        use rand::rngs::StdRng;
        use rand::Rng;
        use rand::SeedableRng;
        let mut rng = StdRng::seed_from_u64(seed);
        (0..len)
            .map(|_| match rng.random_range(0u8..4) {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect()
    }

    fn write_fasta(path: &str, records: &[(&str, &[u8])]) {
        use std::io::Write;
        let mut f = std::fs::File::create(path).unwrap();
        for (name, seq) in records {
            writeln!(f, ">{}", name).unwrap();
            writeln!(f, "{}", std::str::from_utf8(seq).unwrap()).unwrap();
        }
    }

    #[test]
    fn test_decompressor_basic() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        let ref_seq = random_dna(2000, 42);
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let sample_path = dir.path().join("sample.fa");
        let mut sample_seq = ref_seq.clone();
        sample_seq[100] = b'G';
        sample_seq[200] = b'C';
        write_fasta(sample_path.to_str().unwrap(), &[("chr1", &sample_seq)]);

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample("sample1", sample_path.to_str().unwrap())?;
        comp.finish()?;

        // Open with Decompressor.
        let dec = Decompressor::open(&out_path)?;
        assert_eq!(dec.list_samples(), vec!["sample1"]);
        assert!(dec.contains_contig("chr1"));
        assert!(!dec.contains_contig("chr2"));
        Ok(())
    }

    #[test]
    fn test_get_sample_roundtrip() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        let ref_seq = random_dna(2000, 42);
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let sample_path = dir.path().join("sample.fa");
        let mut sample_seq = ref_seq.clone();
        sample_seq[100] = b'G';
        sample_seq[200] = b'C';
        sample_seq[300] = b'T';
        write_fasta(sample_path.to_str().unwrap(), &[("chr1", &sample_seq)]);

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample("sample1", sample_path.to_str().unwrap())?;
        comp.finish()?;

        // Extract sample back.
        let mut dec = Decompressor::open(&out_path)?;
        let mut out_buf = Vec::new();
        dec.get_sample("sample1", &mut out_buf)?;

        let out_str = String::from_utf8(out_buf)?;
        // The output should contain the sample sequence (uppercase, since
        // 2-bit encoding loses case info).
        let expected =
            String::from_utf8(sample_seq.iter().map(|&c| c.to_ascii_uppercase()).collect())
                .unwrap();
        // Check that the sequence appears in the output (after the header line).
        let lines: Vec<&str> = out_str.lines().collect();
        assert!(lines[0].starts_with(">chr1"));
        let seq: String = lines[1..].concat();
        assert_eq!(seq, expected);
        Ok(())
    }

    #[test]
    fn test_get_contig_roundtrip() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        let ref_seq = random_dna(2000, 42);
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let sample_path = dir.path().join("sample.fa");
        let sample_seq = random_dna(2000, 100);
        write_fasta(sample_path.to_str().unwrap(), &[("chr1", &sample_seq)]);

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample("s1", sample_path.to_str().unwrap())?;
        comp.finish()?;

        let mut dec = Decompressor::open(&out_path)?;
        let mut out_buf = Vec::new();
        dec.get_contig("chr1", None, None, "+", &mut out_buf)?;

        let out_str = String::from_utf8(out_buf)?;
        let lines: Vec<&str> = out_str.lines().collect();
        assert!(lines[0].starts_with(">s1"));
        let seq: String = lines[1..].concat();
        let expected =
            String::from_utf8(sample_seq.iter().map(|&c| c.to_ascii_uppercase()).collect())
                .unwrap();
        assert_eq!(seq, expected);
        Ok(())
    }

    #[test]
    fn test_get_contig_slice() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        let ref_seq = random_dna(2000, 42);
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let sample_path = dir.path().join("sample.fa");
        let sample_seq = random_dna(2000, 100);
        write_fasta(sample_path.to_str().unwrap(), &[("chr1", &sample_seq)]);

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample("s1", sample_path.to_str().unwrap())?;
        comp.finish()?;

        let mut dec = Decompressor::open(&out_path)?;
        let mut out_buf = Vec::new();
        dec.get_contig("chr1", Some(10), Some(20), "+", &mut out_buf)?;

        let out_str = String::from_utf8(out_buf)?;
        let lines: Vec<&str> = out_str.lines().collect();
        let seq: String = lines[1..].concat();
        assert_eq!(seq.len(), 10);
        let expected = String::from_utf8(
            sample_seq[10..20]
                .iter()
                .map(|&c| c.to_ascii_uppercase())
                .collect(),
        )
        .unwrap();
        assert_eq!(seq, expected);
        Ok(())
    }

    #[test]
    fn test_get_contig_neg_strand() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        let ref_seq = random_dna(2000, 42);
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let sample_path = dir.path().join("sample.fa");
        let sample_seq = random_dna(2000, 100);
        write_fasta(sample_path.to_str().unwrap(), &[("chr1", &sample_seq)]);

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample("s1", sample_path.to_str().unwrap())?;
        comp.finish()?;

        let mut dec = Decompressor::open(&out_path)?;
        let mut out_buf = Vec::new();
        dec.get_contig("chr1", Some(0), Some(10), "-", &mut out_buf)?;

        let out_str = String::from_utf8(out_buf)?;
        let lines: Vec<&str> = out_str.lines().collect();
        let seq: String = lines[1..].concat();
        let fwd: Vec<u8> = sample_seq[0..10]
            .iter()
            .map(|&c| c.to_ascii_uppercase())
            .collect();
        let expected: Vec<u8> = nt::rev_comp(&fwd).collect();
        assert_eq!(seq.as_bytes(), expected);
        Ok(())
    }

    #[test]
    fn test_sequence_reader_reference_layer() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        let ref_seq = random_dna(2000, 42);
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let out_path = dir.path().join("out.pbit");
        let comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.finish()?;
        let mut dec = Decompressor::open(&out_path)?;
        // Read full reference.
        let seq = dec.read_sequence("chr1", None, None)?;
        assert_eq!(seq.len(), 2000);
        // Read a slice.
        let slice = dec.read_sequence("chr1", Some(10), Some(20))?;
        assert_eq!(slice.len(), 10);
        // Read missing contig.
        assert!(dec.read_sequence("chr2", None, None).is_err());
        Ok(())
    }

    #[test]
    fn test_multi_segment_contig() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        // 5000 bp → 2 segments of 4096 + 904.
        let ref_seq = random_dna(5000, 42);
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let sample_path = dir.path().join("sample.fa");
        let mut sample_seq = ref_seq.clone();
        sample_seq[4500] = b'G';
        write_fasta(sample_path.to_str().unwrap(), &[("chr1", &sample_seq)]);

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample("s1", sample_path.to_str().unwrap())?;
        comp.finish()?;

        let mut dec = Decompressor::open(&out_path)?;
        let mut out_buf = Vec::new();
        dec.get_sample("s1", &mut out_buf)?;
        let out_str = String::from_utf8(out_buf)?;
        let lines: Vec<&str> = out_str.lines().collect();
        let seq: String = lines[1..].concat();
        assert_eq!(seq.len(), 5000);
        let expected =
            String::from_utf8(sample_seq.iter().map(|&c| c.to_ascii_uppercase()).collect())
                .unwrap();
        assert_eq!(seq, expected);
        Ok(())
    }

    #[test]
    fn test_multiple_samples_get_contig() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        let ref_seq = random_dna(1000, 42);
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let s1_path = dir.path().join("s1.fa");
        let s2_path = dir.path().join("s2.fa");
        let s1_seq = random_dna(1000, 100);
        let s2_seq = random_dna(1000, 200);
        write_fasta(s1_path.to_str().unwrap(), &[("chr1", &s1_seq)]);
        write_fasta(s2_path.to_str().unwrap(), &[("chr1", &s2_seq)]);

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample("s1", s1_path.to_str().unwrap())?;
        comp.append_sample("s2", s2_path.to_str().unwrap())?;
        comp.finish()?;

        let mut dec = Decompressor::open(&out_path)?;
        assert_eq!(dec.list_samples().len(), 2);
        let mut out_buf = Vec::new();
        dec.get_contig("chr1", None, None, "+", &mut out_buf)?;
        let out_str = String::from_utf8(out_buf)?;
        // Should have 2 FASTA entries (one per sample).
        let headers: Vec<&str> = out_str.lines().filter(|l| l.starts_with('>')).collect();
        assert_eq!(headers.len(), 2);
        Ok(())
    }
}
