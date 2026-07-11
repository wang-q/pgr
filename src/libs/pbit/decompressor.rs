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

use super::cigar_delta::{apply_cigar, unpack_cigar};
use super::collection::{Collection, SegmentDesc};
use super::format::{
    read_ref_index, read_u32_le, DeltaEncoding, DeltaMeta, PbitFooter, PbitHeader, RefGroupEntry,
};
use super::segment::Segment;

/// Decompressor for a `.pbit` archive.
pub struct Decompressor<R: Read + Seek> {
    reader: R,
    header: PbitHeader,
    footer: PbitFooter,
    ref_groups: Vec<RefGroupEntry>,
    /// contig name → Vec<ref_group_id> (reference segments, ordered).
    contig_groups: IndexMap<String, Vec<u32>>,
    /// All contig names appearing in any sample's collection (for
    /// `contains_contig`).
    contig_set: HashSet<String>,
    collection: Collection,
    /// delta_meta[ref_group_id][delta_id] → header info (no packed data).
    /// Used by `get_contig` to compute segment coordinates for smart slice
    /// selection (skip non-overlapping segments).
    delta_meta: Vec<Vec<DeltaMeta>>,
    /// delta_offsets[ref_group_id][delta_id] → file offset of the delta's
    /// 10-byte header (followed by `packed_size` bytes).
    delta_offsets: Vec<Vec<u64>>,
    /// LRU cache: ref_group_id → decoded reference segment DNA (ASCII).
    ref_cache: LruCache<u32, Vec<u8>>,
    /// LRU cache: (ref_group_id, delta_id, ref_start, ref_end) → decoded raw
    /// sample segment. ref_start/ref_end are part of the key because
    /// CIGAR-mode deltas decode against a ref slice [ref_start, ref_end); two
    /// segments sharing a delta_id (via packed_data dedup) but with different
    /// ref slices produce different outputs.
    delta_cache: LruCache<(u32, u32, u32, u32), Vec<u8>>,
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

impl Decompressor<std::io::Cursor<Vec<u8>>> {
    /// Open and read entire file into memory (mirrors `TwoBitFile::open_and_read`).
    pub fn open_and_read<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mut buf = Vec::new();
        std::fs::File::open(&path)
            .with_context(|| format!("failed to open pbit file: {}", path.as_ref().display()))?
            .read_to_end(&mut buf)?;
        Self::new(std::io::Cursor::new(buf))
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

        // Scan delta data: read each delta's 10-byte header, build delta_meta
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

        // Read sample index (collection, flate2-compressed). The collection
        // spans [sample_index_offset, footer_start) where footer_start =
        // file_size - 24.
        let file_size = reader.seek(SeekFrom::End(0))?;
        let footer_start = file_size
            .checked_sub(24)
            .ok_or_else(|| anyhow!("pbit file too small: {} bytes", file_size))?;
        if footer.sample_index_offset > footer_start {
            return Err(anyhow!(
                "sample_index_offset {} exceeds footer start {}",
                footer.sample_index_offset,
                footer_start
            ));
        }
        let collection_len = footer_start - footer.sample_index_offset;
        reader.seek(SeekFrom::Start(footer.sample_index_offset))?;
        let mut compressed = vec![0u8; collection_len as usize];
        reader.read_exact(&mut compressed)?;
        let collection = Collection::deserialize(&compressed)?;

        // Validate sample_count consistency.
        if header.sample_count != collection.samples.len() as u32 {
            return Err(anyhow!(
                "sample_count mismatch: header={}, collection={}",
                header.sample_count,
                collection.samples.len()
            ));
        }

        // Build contig_set from collection.
        let mut contig_set: HashSet<String> = HashSet::new();
        for (_sample, contigs) in &collection.samples {
            for cs in contigs {
                contig_set.insert(cs.contig_name.clone());
            }
        }

        let min_match_len = header.min_match_len;

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

    /// Return the footer (for `Compressor::open_for_append`).
    pub fn footer(&self) -> &PbitFooter {
        &self.footer
    }

    /// Return an owned clone of the collection (for `Compressor::open_for_append`).
    pub fn collection_clone(&self) -> Collection {
        self.collection.clone()
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

    /// Read a delta's packed data and decode it (LZ-diff or CIGAR depending
    /// on `encoding`) against the reference segment. Uses an LRU cache.
    /// CIGAR-encoded deltas use `seg.ref_start` / `seg.ref_end` to slice the
    /// reference; LZ-diff deltas ignore them. Both encodings store
    /// gzip-compressed packed_data, but the decompression path differs:
    /// LZ-diff uses `flate2::read::GzDecoder` then `Segment::get`; CIGAR uses
    /// `unpack_cigar` (which includes its own gzip decompression).
    fn decode_delta(&mut self, seg: &SegmentDesc) -> Result<Vec<u8>> {
        let key = (seg.ref_group_id, seg.delta_id, seg.ref_start, seg.ref_end);
        if let Some(cached) = self.delta_cache.get(&key) {
            return Ok(cached.clone());
        }

        // Validate SegmentDesc indices against the in-file metadata. These
        // values come from the (potentially corrupted) sample index, so they
        // must be bounds-checked before use to avoid panics.
        let gid = seg.ref_group_id as usize;
        let did = seg.delta_id as usize;
        if gid >= self.ref_groups.len() {
            anyhow::bail!(
                "decode_delta: ref_group_id {} out of range ({})",
                seg.ref_group_id,
                self.ref_groups.len()
            );
        }
        if did >= self.delta_offsets[gid].len() {
            anyhow::bail!(
                "decode_delta: delta_id {} out of range ({}) for ref_group {}",
                seg.delta_id,
                self.delta_offsets[gid].len(),
                seg.ref_group_id
            );
        }

        // Read reference segment.
        let ref_dna = self.read_ref_segment(seg.ref_group_id)?;

        // Read packed delta data. The 10-byte header was already scanned at
        // construction and cached in self.delta_meta, so seek past it.
        let offset = self.delta_offsets[gid][did];
        let meta = self.delta_meta[gid][did];
        self.reader.seek(SeekFrom::Start(offset + 10))?;
        let mut packed = vec![0u8; meta.packed_size as usize];
        self.reader.read_exact(&mut packed)?;

        // Decode by encoding type.
        let decoded = match meta.encoding {
            DeltaEncoding::LzDiff => {
                // LZ-diff: packed_data is flate2-compressed raw delta.
                let mut decoder = flate2::read::GzDecoder::new(&packed[..]);
                let mut delta = Vec::new();
                decoder.read_to_end(&mut delta)?;
                let mut lz = Segment::new(self.min_match_len);
                lz.prepare(&ref_dna);
                lz.get(&delta)?
            }
            DeltaEncoding::Cigar => {
                // CIGAR: packed_data is pack_cigar output (includes its own gzip).
                if seg.ref_start >= seg.ref_end {
                    anyhow::bail!(
                        "decode_delta: invalid CIGAR reference interval [{}; {})",
                        seg.ref_start,
                        seg.ref_end
                    );
                }
                if (seg.ref_end as usize) > ref_dna.len() {
                    anyhow::bail!(
                        "decode_delta: ref_end {} > ref segment length {}",
                        seg.ref_end,
                        ref_dna.len()
                    );
                }
                let (ops, xi_bases) = unpack_cigar(&packed)?;
                let ref_slice = &ref_dna[seg.ref_start as usize..seg.ref_end as usize];
                apply_cigar(ref_slice, &ops, &xi_bases)?
            }
        };

        // Apply reverse-complement if needed.
        let final_decoded = if meta.is_rev_comp {
            nt::rev_comp(&decoded).collect()
        } else {
            decoded
        };

        self.delta_cache.put(key, final_decoded.clone());
        Ok(final_decoded)
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
            // Extract raw_lengths first to release the immutable borrow on
            // self.delta_meta before calling self.decode_delta (mutable).
            let seg_lens: Vec<usize> = segments
                .iter()
                .map(|seg| {
                    let gid = seg.ref_group_id as usize;
                    let did = seg.delta_id as usize;
                    let meta_row = self.delta_meta.get(gid).ok_or_else(|| {
                        anyhow::anyhow!(
                            "get_contig: ref_group_id {} out of range ({})",
                            seg.ref_group_id,
                            self.delta_meta.len()
                        )
                    })?;
                    let meta = meta_row.get(did).ok_or_else(|| {
                        anyhow::anyhow!(
                            "get_contig: delta_id {} out of range ({}) for ref_group {}",
                            seg.delta_id,
                            meta_row.len(),
                            seg.ref_group_id
                        )
                    })?;
                    Ok::<usize, anyhow::Error>(meta.raw_length as usize)
                })
                .collect::<Result<Vec<_>>>()?;
            let total_len: usize = seg_lens.iter().sum();

            // Clamp [s, e) to [0, total_len].
            let s = start.unwrap_or(0).min(total_len);
            let e = end.unwrap_or(total_len).min(total_len);
            if s >= e {
                continue;
            }

            // Decode only segments overlapping [s, e) (smart selection, like
            // the reference layer's read_sequence).
            let mut result = Vec::new();
            let mut offset: usize = 0;
            for (seg, &seg_len) in segments.iter().zip(seg_lens.iter()) {
                let seg_end = offset + seg_len;
                if seg_end > s && offset < e {
                    let decoded = self.decode_delta(seg)?;
                    anyhow::ensure!(
                        decoded.len() == seg_len,
                        "decoded segment length {} does not match metadata raw_length {} \
                         for sample '{}' contig '{}' (archive may be corrupt)",
                        decoded.len(),
                        seg_len,
                        sample,
                        contig
                    );
                    let local_start = s.saturating_sub(offset).min(decoded.len());
                    let local_end = (e - offset).min(seg_len).min(decoded.len());
                    if local_start < local_end {
                        result.extend_from_slice(&decoded[local_start..local_end]);
                    }
                }
                offset = seg_end;
                if offset >= e {
                    break;
                }
            }

            // Apply reverse-complement if needed.
            let seq_bytes: Vec<u8> = if strand == "-" {
                nt::rev_comp(&result).collect()
            } else {
                result
            };

            // Write FASTA header.
            let header = match (start, end) {
                (Some(_), Some(_)) => {
                    format!(">{} {}:{}-{}({})", sample, contig, s + 1, e, strand)
                }
                _ if strand == "-" => format!(">{} {}(-)", sample, contig),
                _ => format!(">{} {}", sample, contig),
            };
            writeln!(out, "{}", header)?;
            write_fasta_seq(out, &seq_bytes, line_width)?;
        }
        Ok(())
    }

    /// Extract all contigs of a single sample, writing FASTA entries.
    pub fn get_sample(&mut self, sample: &str, out: &mut impl Write) -> Result<()> {
        let line_width = 60;

        // Collect (contig_name, segments) pairs first to release the immutable
        // borrow on self.collection before calling self.decode_delta().
        let contig_segs: Vec<(String, Vec<SegmentDesc>)> = match self.collection.samples.get(sample)
        {
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
                let decoded = self.decode_delta(seg)?;

                // Validate decoded length against cached metadata (indices come
                // from the potentially corrupted sample index).
                let gid = seg.ref_group_id as usize;
                let did = seg.delta_id as usize;
                let expected = self
                    .delta_meta
                    .get(gid)
                    .and_then(|row| row.get(did))
                    .map(|m| m.raw_length as usize)
                    .ok_or_else(|| {
                        anyhow!(
                            "get_sample: ref_group_id {} or delta_id {} out of range",
                            seg.ref_group_id,
                            seg.delta_id
                        )
                    })?;
                anyhow::ensure!(
                    decoded.len() == expected,
                    "decoded segment length {} does not match metadata raw_length {} \
                     for sample '{}' contig '{}' (archive may be corrupt)",
                    decoded.len(),
                    expected,
                    sample,
                    contig_name
                );

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

    #[test]
    fn test_get_contig_corrupt_raw_length_returns_error() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        let ref_seq = random_dna(1000, 42);
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let sample_path = dir.path().join("sample.fa");
        let sample_seq = random_dna(1000, 100);
        write_fasta(sample_path.to_str().unwrap(), &[("chr1", &sample_seq)]);

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample("s1", sample_path.to_str().unwrap())?;
        comp.finish()?;

        // Patch the first delta's raw_length to be larger than the decoded segment.
        // Delta data layout: delta_data_offset + 4 (ref_group_count) + 4 (delta_count)
        // + 1 (is_rev_comp) -> raw_length u32.
        let mut file = std::fs::File::open(&out_path)?;
        file.seek(SeekFrom::End(-24))?;
        let mut footer_buf = [0u8; 24];
        file.read_exact(&mut footer_buf)?;
        let delta_data_offset = u64::from_le_bytes([
            footer_buf[8],
            footer_buf[9],
            footer_buf[10],
            footer_buf[11],
            footer_buf[12],
            footer_buf[13],
            footer_buf[14],
            footer_buf[15],
        ]);
        let raw_length_offset = delta_data_offset + 4 + 4 + 1;
        let mut file = std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .open(&out_path)?;
        file.seek(SeekFrom::Start(raw_length_offset))?;
        // Write a deliberately wrong raw_length (2000 instead of 1000).
        file.write_all(&2000u32.to_le_bytes())?;
        drop(file);

        let mut dec = Decompressor::open(&out_path)?;
        let mut out_buf = Vec::new();
        let res = dec.get_contig("chr1", None, None, "+", &mut out_buf);
        assert!(res.is_err());
        let err = res.unwrap_err().to_string();
        assert!(
            err.contains("does not match metadata raw_length"),
            "unexpected error: {}",
            err
        );
        Ok(())
    }

    #[test]
    fn test_get_sample_corrupt_raw_length_returns_error() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        let ref_seq = random_dna(1000, 42);
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let sample_path = dir.path().join("sample.fa");
        let sample_seq = random_dna(1000, 100);
        write_fasta(sample_path.to_str().unwrap(), &[("chr1", &sample_seq)]);

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample("s1", sample_path.to_str().unwrap())?;
        comp.finish()?;

        // Patch the first delta's raw_length to be larger than the decoded segment.
        let mut file = std::fs::File::open(&out_path)?;
        file.seek(SeekFrom::End(-24))?;
        let mut footer_buf = [0u8; 24];
        file.read_exact(&mut footer_buf)?;
        let delta_data_offset = u64::from_le_bytes([
            footer_buf[8],
            footer_buf[9],
            footer_buf[10],
            footer_buf[11],
            footer_buf[12],
            footer_buf[13],
            footer_buf[14],
            footer_buf[15],
        ]);
        let raw_length_offset = delta_data_offset + 4 + 4 + 1;
        let mut file = std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .open(&out_path)?;
        file.seek(SeekFrom::Start(raw_length_offset))?;
        // Write a deliberately wrong raw_length (2000 instead of 1000).
        file.write_all(&2000u32.to_le_bytes())?;
        drop(file);

        let mut dec = Decompressor::open(&out_path)?;
        let mut out_buf = Vec::new();
        let res = dec.get_sample("s1", &mut out_buf);
        assert!(res.is_err());
        let err = res.unwrap_err().to_string();
        assert!(
            err.contains("does not match metadata raw_length"),
            "unexpected error: {}",
            err
        );
        Ok(())
    }
}
