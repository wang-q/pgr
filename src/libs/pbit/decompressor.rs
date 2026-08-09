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

/// (contig_name, segments, mask_blocks) gathered for one sample's contigs.
type SampleContig = (String, Vec<SegmentDesc>, Vec<(u32, u32)>);

/// Parse the flate2-compressed PAF recovery block written by
/// `Compressor::write_paf_data` (v1009).
fn read_paf_data(data: &[u8]) -> Result<Vec<PafRecovery>> {
    let raw = crate::libs::bgzf::gzip_decompress(data, MAX_DELTA_UNCOMPRESSED)?;
    if raw.len() > MAX_DELTA_UNCOMPRESSED {
        anyhow::bail!(
            "PAF recovery data exceeds maximum uncompressed size {} bytes",
            MAX_DELTA_UNCOMPRESSED
        );
    }
    let mut cursor = std::io::Cursor::new(raw);
    let sample_count = read_u32_le(&mut cursor)? as usize;
    let mut out = Vec::with_capacity(sample_count.min(1024));
    for _ in 0..sample_count {
        let sample = read_string(&mut cursor)?;
        let big_count = read_u32_le(&mut cursor)? as usize;
        let mut big_ms = Vec::with_capacity(big_count.min(1 << 16));
        for _ in 0..big_count {
            let record_id = read_u32_le(&mut cursor)?;
            let ms = read_u32_le(&mut cursor)? as i32;
            big_ms.push((record_id, ms));
        }
        let small_count = read_u32_le(&mut cursor)? as usize;
        let mut small = Vec::with_capacity(small_count.min(1 << 16));
        for _ in 0..small_count {
            small.push(read_string(&mut cursor)?);
        }
        out.push((sample, big_ms, small));
    }
    Ok(out)
}
use super::format::{
    read_ref_index, read_ref_table, read_string, read_u32_le, DeltaEncoding, DeltaMeta,
    PafRecovery, PbitFooter, PbitHeader, RefGroupEntry, RefTableEntry, MAX_DELTAS_PER_GROUP,
    MAX_DELTA_UNCOMPRESSED, MAX_PACKED_SIZE, MAX_REF_GROUPS,
};
use super::segment::Segment;

/// FASTA line wrap width for output.
const FASTA_LINE_WIDTH: usize = 60;

/// Decompressor for a `.pbit` archive.
pub struct Decompressor<R: Read + Seek> {
    reader: R,
    header: PbitHeader,
    footer: PbitFooter,
    ref_groups: Vec<RefGroupEntry>,
    /// Per-reference metadata (name, group range, embedded-index offsets).
    ref_meta: Vec<RefTableEntry>,
    /// contig name → Vec<ref_group_id> (reference segments, ordered).
    contig_groups: IndexMap<String, Vec<u32>>,
    /// ref_group_id → global reference offset of the segment's first base
    /// (v1007 CIGAR deltas reference reference-file-global coordinates).
    seg_starts: Vec<u64>,
    /// Total reference length across all reference files.
    ref_total_len: u64,
    /// All contig names appearing in any sample's collection (for
    /// `contains_contig`).
    contig_set: HashSet<String>,
    collection: Collection,
    /// PAF recovery data per sample (v1009): (sample, big-chain ms table,
    /// verbatim small-chain PAF rows).
    paf_data: Vec<PafRecovery>,
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

        // Guard against a malicious/corrupt archive: `min_match_len` drives
        // `LzDiff::prepare`'s `reference.resize(len + key_len)` padding, so a
        // value far larger than segment_size would trigger a multi-GB allocation
        // per reference segment during extraction. A match cannot span more than
        // one segment, so min_match_len > segment_size is always invalid.
        if header.min_match_len > header.segment_size {
            return Err(anyhow!(
                "archive min_match_len {} exceeds segment_size {}; corrupt or malicious",
                header.min_match_len,
                header.segment_size
            ));
        }
        // `min_match_len` drives `LzDiff::prepare`'s `reference.resize(len + key_len)`
        // padding (key_len ≈ min_match_len), so a huge value would force a
        // multi-GB allocation per decoded segment even when `segment_size` is
        // also huge enough to pass the segment_size check above. A match cannot
        // span more than one segment, so cap it at the per-segment payload bound.
        if header.min_match_len as usize > MAX_PACKED_SIZE {
            return Err(anyhow!(
                "archive min_match_len {} exceeds maximum {}; corrupt or malicious",
                header.min_match_len,
                MAX_PACKED_SIZE
            ));
        }
        // A zero segment_size / kmer_len is invalid for any well-formed archive
        // (create validates both as positive). `open_for_append` reuses these
        // values to re-segment sample/reference FASTAs, and `segment_sequence`
        // calls `chunks(0)` / `detect_rev_comp` calls `windows(0)`, both of
        // which panic. Reject here so a crafted archive cannot crash `append` /
        // `append-ref` (Zero Panic).
        if header.segment_size == 0 {
            return Err(anyhow!(
                "archive segment_size must be positive; corrupt or malicious"
            ));
        }
        if header.kmer_len == 0 {
            return Err(anyhow!(
                "archive kmer_len must be positive; corrupt or malicious"
            ));
        }
        // Bound ref_group_count: it drives `Vec::with_capacity` for the delta
        // scan below, so a malicious value would force a multi-GB allocation
        // before any data is validated.
        if header.ref_group_count as usize > MAX_REF_GROUPS {
            return Err(anyhow!(
                "archive ref_group_count {} exceeds maximum {}; corrupt or malicious",
                header.ref_group_count,
                MAX_REF_GROUPS
            ));
        }

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
        let ref_meta = read_ref_table(&mut reader)?;
        if ref_meta.is_empty() {
            return Err(anyhow!("reference table is empty"));
        }

        // Build contig_groups (contig name → ref_group_ids in order).
        let mut contig_groups: IndexMap<String, Vec<u32>> = IndexMap::new();
        for (i, entry) in ref_groups.iter().enumerate() {
            contig_groups
                .entry(entry.contig_name.clone())
                .or_default()
                .push(i as u32);
        }

        // Build the reference coordinate index (global offset per segment).
        let mut seg_starts = Vec::with_capacity(ref_groups.len());
        let mut ref_pos = 0u64;
        for entry in &ref_groups {
            seg_starts.push(ref_pos);
            reader.seek(SeekFrom::Start(entry.segment_offset))?;
            let dna_size = read_u32_le(&mut reader)? as u64;
            ref_pos += dna_size;
        }
        let ref_total_len = ref_pos;

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
            if delta_count > MAX_DELTAS_PER_GROUP {
                return Err(anyhow!(
                    "delta_count {} exceeds maximum {}; corrupt or malicious",
                    delta_count,
                    MAX_DELTAS_PER_GROUP
                ));
            }
            let mut metas = Vec::with_capacity(delta_count);
            let mut offsets = Vec::with_capacity(delta_count);
            for _ in 0..delta_count {
                let offset = reader.stream_position()?;
                let meta = DeltaMeta::read_header(&mut reader)?;
                // An inflated packed_size would trigger a multi-GB allocation
                // in `decode_delta`; reject it while scanning (before any data
                // is read into memory).
                if meta.packed_size as usize > MAX_PACKED_SIZE {
                    return Err(anyhow!(
                        "delta packed_size {} exceeds maximum {}; corrupt or malicious",
                        meta.packed_size,
                        MAX_PACKED_SIZE
                    ));
                }
                metas.push(meta);
                offsets.push(offset);
                // Skip the packed data.
                reader.seek(SeekFrom::Current(meta.packed_size as i64))?;
            }
            delta_meta.push(metas);
            delta_offsets.push(offsets);
        }

        // Read sample index (collection, flate2-compressed) and the PAF
        // recovery data. The collection spans [sample_index_offset,
        // paf_data_offset) and PAF data spans [paf_data_offset, footer_start)
        // where footer_start = file_size - 32 (v1009 footer).
        let file_size = reader.seek(SeekFrom::End(0))?;
        let footer_start = file_size
            .checked_sub(32)
            .ok_or_else(|| anyhow!("pbit file too small: {} bytes", file_size))?;
        if footer.sample_index_offset > footer_start {
            return Err(anyhow!(
                "sample_index_offset {} exceeds footer start {}",
                footer.sample_index_offset,
                footer_start
            ));
        }
        let collection_end = footer.paf_data_offset.min(footer_start);
        let collection_len = collection_end - footer.sample_index_offset;
        reader.seek(SeekFrom::Start(footer.sample_index_offset))?;
        let mut compressed = vec![0u8; collection_len as usize];
        reader.read_exact(&mut compressed)?;
        let collection = Collection::deserialize(&compressed)?;

        // PAF recovery data (v1009; empty for archives without it).
        let paf_data = if footer.paf_data_offset < footer_start {
            let paf_len = footer_start - footer.paf_data_offset;
            reader.seek(SeekFrom::Start(footer.paf_data_offset))?;
            let mut paf_compressed = vec![0u8; paf_len as usize];
            reader.read_exact(&mut paf_compressed)?;
            read_paf_data(&paf_compressed)?
        } else {
            Vec::new()
        };

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
            ref_meta,
            contig_groups,
            seg_starts,
            ref_total_len,
            contig_set,
            collection,
            paf_data,
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

    /// Global reference offset of each reference segment's first base.
    pub fn ref_seg_starts(&self) -> &[u64] {
        &self.seg_starts
    }

    /// Return the reference contig name, the contig-relative start of the
    /// given reference segment, and the contig's total length. Used by
    /// `pbit to-paf` to project global reference coordinates back to contig
    /// coordinates.
    pub fn ref_group_location(&self, ref_group_id: u32) -> Option<(String, u32, u32)> {
        let entry = self.ref_groups.get(ref_group_id as usize)?;
        let segs = self.contig_groups.get(&entry.contig_name)?;
        let idx = segs.iter().position(|&g| g == ref_group_id)?;
        let seg_size = self.header.segment_size;
        let start = idx as u32 * seg_size;
        // Contig total length = all full segments plus the actual length of
        // the last (possibly shorter) segment.
        let last_gid = segs[segs.len() - 1];
        let g_start = self.seg_starts[last_gid as usize];
        let g_end = if last_gid as usize + 1 < self.seg_starts.len() {
            self.seg_starts[last_gid as usize + 1]
        } else {
            self.ref_total_len
        };
        let total = (segs.len() - 1) as u32 * seg_size + (g_end - g_start) as u32;
        Some((entry.contig_name.clone(), start, total))
    }

    /// Read one segment's delta metadata and packed payload (without
    /// decoding). Used by `pbit to-paf` to recover embedded CIGAR alignments.
    pub fn segment_payload(&mut self, seg: &SegmentDesc) -> Result<(DeltaMeta, Vec<u8>)> {
        let gid = seg.ref_group_id as usize;
        let did = seg.delta_id as usize;
        let meta = self
            .delta_meta
            .get(gid)
            .and_then(|row| row.get(did))
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "segment_payload: ref_group {} / delta {} out of range",
                    seg.ref_group_id,
                    seg.delta_id
                )
            })?;
        let offset = self.delta_offsets[gid][did];
        self.reader.seek(SeekFrom::Start(offset + 10))?;
        let mut packed = vec![0u8; meta.packed_size as usize];
        self.reader.read_exact(&mut packed)?;
        Ok((*meta, packed))
    }

    /// Per-reference metadata (name, group range, embedded-index offsets).
    pub fn ref_table(&self) -> &[RefTableEntry] {
        &self.ref_meta
    }

    /// Return the header (for `stat` overview).
    pub fn header(&self) -> &PbitHeader {
        &self.header
    }

    /// Return the collection (for `stat --contigs`).
    pub fn collection(&self) -> &Collection {
        &self.collection
    }

    /// PAF recovery data per sample (v1009): big-chain (record_id, ms) tables
    /// and verbatim small-chain PAF rows.
    pub fn paf_data(&self) -> &[PafRecovery] {
        &self.paf_data
    }

    /// Count referenced deltas by encoding type (LzDiff/Cigar/Raw/Identity).
    pub fn delta_encoding_counts(&self) -> [usize; 4] {
        let mut counts = [0usize; 4];
        for contigs in self.collection.samples.values() {
            for contig in contigs {
                for seg in &contig.segments {
                    if let Some(row) = self.delta_meta.get(seg.ref_group_id as usize) {
                        if let Some(meta) = row.get(seg.delta_id as usize) {
                            counts[meta.encoding as usize] += 1;
                        }
                    }
                }
            }
        }
        counts
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

    /// Read an arbitrary reference interval `[start, end)` in
    /// reference-file-global coordinates (v1007 CIGAR deltas). May span
    /// multiple reference segments; the interval is stitched from segment
    /// slices.
    pub fn read_ref_interval(&mut self, start: u64, end: u64) -> Result<Vec<u8>> {
        if start >= end {
            return Ok(Vec::new());
        }
        if end > self.ref_total_len {
            anyhow::bail!(
                "read_ref_interval: end {} exceeds total reference length {}",
                end,
                self.ref_total_len
            );
        }
        let mut idx = match self.seg_starts.binary_search(&start) {
            Ok(i) => i,
            Err(0) => {
                anyhow::bail!("read_ref_interval: start {} before first segment", start)
            }
            Err(i) => i - 1,
        };
        let mut out = Vec::with_capacity((end - start) as usize);
        let mut cur = start;
        while cur < end {
            let seg_start = self.seg_starts[idx];
            let seg_end = if idx + 1 < self.seg_starts.len() {
                self.seg_starts[idx + 1]
            } else {
                self.ref_total_len
            };
            let seg = self.read_ref_segment(idx as u32)?;
            let lo = (cur - seg_start) as usize;
            let hi = ((end - seg_start) as usize).min((seg_end - seg_start) as usize);
            if lo > seg.len() || hi > seg.len() {
                anyhow::bail!(
                    "read_ref_interval: segment {idx} slice [{lo}, {hi}) out of range (len {})",
                    seg.len()
                );
            }
            out.extend_from_slice(&seg[lo..hi]);
            cur = seg_start + hi as u64;
            idx += 1;
        }
        Ok(out)
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
                // Read the owning reference segment (LZ-diff is segment-relative).
                let ref_dna = self.read_ref_segment(seg.ref_group_id)?;
                // LZ-diff: packed_data is flate2-compressed raw delta. Bound the
                // decompressed size to reject gzip bombs (an attacker could
                // otherwise expand a tiny payload into a multi-GB allocation).
                let delta = crate::libs::bgzf::gzip_decompress(&packed, MAX_DELTA_UNCOMPRESSED)?;
                if delta.len() > MAX_DELTA_UNCOMPRESSED {
                    anyhow::bail!(
                        "LZ-diff delta decompressed size exceeds maximum {} bytes",
                        MAX_DELTA_UNCOMPRESSED
                    );
                }
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
                let (ops, xi_bases) = unpack_cigar(&packed)?;
                // ref_start/ref_end are reference-file-global coordinates.
                let ref_dna = self.read_ref_interval(seg.ref_start as u64, seg.ref_end as u64)?;
                apply_cigar(&ref_dna, &ops, &xi_bases)?
            }
            DeltaEncoding::Raw => {
                // Verbatim segment as a standard 2bit record; the reference
                // segment is not used.
                let mut cursor = std::io::Cursor::new(&packed[..]);
                let seq = read_2bit_record(&mut cursor, false, None, None, true)?;
                seq.into_bytes()
            }
            DeltaEncoding::Identity => {
                // Zero-payload pointer: the sample segment is exactly the
                // reference interval (rev-comp applied below).
                if seg.ref_start >= seg.ref_end {
                    anyhow::bail!(
                        "decode_delta: invalid Identity reference interval [{}; {})",
                        seg.ref_start,
                        seg.ref_end
                    );
                }
                if seg.ref_end - seg.ref_start != meta.raw_length {
                    anyhow::bail!(
                        "decode_delta: Identity interval length {} does not match raw_length {}",
                        seg.ref_end - seg.ref_start,
                        meta.raw_length
                    );
                }
                self.read_ref_interval(seg.ref_start as u64, seg.ref_end as u64)?
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
    ///
    /// Returns the number of FASTA entries written (0 means the contig exists
    /// but the requested `[start, end)` slice is empty after clamping — e.g.
    /// a range entirely beyond the sequence length). Callers can use this to
    /// warn about silently-empty output.
    pub fn get_contig(
        &mut self,
        contig: &str,
        start: Option<usize>,
        end: Option<usize>,
        strand: &str,
        out: &mut impl Write,
    ) -> Result<usize> {
        let line_width = FASTA_LINE_WIDTH;
        let mut written = 0usize;

        // Collect (sample_name, segments, mask_blocks) first to release the
        // immutable borrow on self.collection before calling self.decode_delta().
        // We must clone the segments (and mask_blocks) (not just borrow) to
        // fully release the immutable borrow on self.collection.
        let sample_segs: Vec<SampleContig> = self
            .collection
            .samples
            .iter()
            .filter_map(|(s, contigs)| {
                contigs
                    .iter()
                    .find(|c| c.contig_name == contig)
                    .map(|cs| (s.clone(), cs.segments.clone(), cs.mask_blocks.clone()))
            })
            .collect();

        for (sample, segments, mask_blocks) in sample_segs {
            // Extract raw_lengths first to release the immutable borrow on
            // self.delta_meta before calling self.decode_delta (mutable).
            // v1007 segments carry explicit contig offsets; order by q_start
            // and compute each segment's contig interval [q_start, q_start+len).
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
            let mut order: Vec<usize> = (0..segments.len()).collect();
            order.sort_by_key(|&i| segments[i].q_start);
            let total_len: usize = order
                .iter()
                .map(|&i| {
                    let seg = &segments[i];
                    (seg.q_start as usize).saturating_add(seg_lens[i])
                })
                .max()
                .unwrap_or(0);

            // Clamp [s, e) to [0, total_len].
            let s = start.unwrap_or(0).min(total_len);
            let e = end.unwrap_or(total_len).min(total_len);
            if s >= e {
                continue;
            }

            // Decode only segments overlapping [s, e) (smart selection, like
            // the reference layer's read_sequence).
            let mut result = Vec::new();
            for &i in &order {
                let seg = &segments[i];
                let seg_start = seg.q_start as usize;
                let seg_end = seg_start + seg_lens[i];
                if seg_end > s && seg_start < e {
                    let decoded = self.decode_delta(seg)?;
                    anyhow::ensure!(
                        decoded.len() == seg_lens[i],
                        "decoded segment length {} does not match metadata raw_length {} \
                         for sample '{}' contig '{}' (archive may be corrupt)",
                        decoded.len(),
                        seg_lens[i],
                        sample,
                        contig
                    );
                    let local_start = s.saturating_sub(seg_start).min(decoded.len());
                    let local_end = (e - seg_start).min(seg_lens[i]).min(decoded.len());
                    if local_start < local_end {
                        result.extend_from_slice(&decoded[local_start..local_end]);
                    }
                }
                if seg_start >= e {
                    break;
                }
            }

            // Apply soft-mask (lowercase) intervals so `some`/`range`
            // extraction is lossless, consistent with `get_sample` (v1005).
            // mask_blocks use 0-based forward contig coords; map to the
            // slice-local [s, e) via offset `s`. Do this on the forward slice
            // before any reverse-complement (rev_comp preserves case).
            apply_mask_blocks_at(&mut result, &mask_blocks, s);

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
            written += 1;
        }
        Ok(written)
    }

    /// Extract all contigs of a single sample, writing FASTA entries.
    pub fn get_sample(&mut self, sample: &str, out: &mut impl Write) -> Result<()> {
        let line_width = FASTA_LINE_WIDTH;

        // Collect (contig_name, segments, mask_blocks) first to release the
        // immutable borrow on self.collection before calling self.decode_delta().
        let contig_segs: Vec<SampleContig> = match self.collection.samples.get(sample) {
            Some(c) => c
                .iter()
                .map(|cs| {
                    (
                        cs.contig_name.clone(),
                        cs.segments.clone(),
                        cs.mask_blocks.clone(),
                    )
                })
                .collect(),
            None => {
                return Err(anyhow!("sample '{}' not found in archive", sample));
            }
        };

        for (contig_name, segments, mask_blocks) in contig_segs {
            let mut full_seq = Vec::new();
            // Segments may be non-contiguous in insertion order (v1007 mixed
            // CIGAR/Raw encoding); order by sample-contig offset and require a
            // gapless tiling so reconstruction is exact.
            let mut ordered = segments.clone();
            ordered.sort_by_key(|s| s.q_start);
            let mut cur = 0u32;
            for seg in &ordered {
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

                if seg.q_start != cur {
                    anyhow::bail!(
                        "gap in sample '{}' contig '{}': segment at {} expected {} (archive may be corrupt)",
                        sample,
                        contig_name,
                        seg.q_start,
                        cur
                    );
                }
                cur = seg.q_start + decoded.len() as u32;
                full_seq.extend_from_slice(&decoded);
            }
            // Restore soft-mask (lowercase) intervals so reconstruction is
            // lossless (v1005; mask_blocks use contig-level 0-based coords).
            apply_mask_blocks(&mut full_seq, &mask_blocks);
            writeln!(out, ">{}", contig_name)?;
            write_fasta_seq(out, &full_seq, line_width)?;
        }
        Ok(())
    }
}

/// Apply soft-mask intervals (lowercase) to a reconstructed sample sequence.
fn apply_mask_blocks(seq: &mut [u8], mask_blocks: &[(u32, u32)]) {
    apply_mask_blocks_at(seq, mask_blocks, 0);
}

/// Apply soft-mask intervals (lowercase) to a slice of a contig sequence.
/// `slice_start` is the 0-based forward offset of `seq` within its contig;
/// mask blocks use 0-based forward contig coordinates.
fn apply_mask_blocks_at(seq: &mut [u8], mask_blocks: &[(u32, u32)], slice_start: usize) {
    for &(start, size) in mask_blocks {
        let b = start as usize;
        let e = b + size as usize;
        let lo = b.saturating_sub(slice_start).min(seq.len());
        let hi = e.saturating_sub(slice_start).min(seq.len());
        if lo < hi {
            for c in &mut seq[lo..hi] {
                *c = c.to_ascii_lowercase();
            }
        }
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
    fn test_decompressor_rejects_huge_ref_group_count() -> Result<()> {
        // A malicious archive with a tiny body but an absurd ref_group_count
        // must be rejected up front (valid magic/version), not crash on a
        // multi-GB allocation.
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("bad.pbit");
        let header = PbitHeader::new(4096, 15, 18, MAX_REF_GROUPS as u32 + 1, 0);
        let footer = PbitFooter {
            ref_index_offset: 36,
            delta_data_offset: 36,
            sample_index_offset: 36,
            paf_data_offset: 36,
        };
        let mut buf = Vec::new();
        header.write_to(&mut buf)?;
        footer.write_to(&mut buf)?;
        std::fs::write(&path, &buf)?;

        let res = Decompressor::open(&path);
        assert!(res.is_err(), "expected rejection for huge ref_group_count");
        let err = res.err().unwrap().to_string();
        assert!(err.contains("ref_group_count"), "unexpected error: {}", err);
        Ok(())
    }

    #[test]
    fn test_decompressor_rejects_zero_segment_size() -> Result<()> {
        // A malicious archive with segment_size=0 must be rejected up front:
        // `open_for_append` re-segments FASTAs with this value and
        // `segment_sequence` calls `chunks(0)`, which panics (Zero Panic).
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("bad.pbit");
        // min_match_len=0 so the earlier `min_match_len > segment_size` and
        // `> MAX_PACKED_SIZE` checks pass, isolating the segment_size==0 check.
        let header = PbitHeader::new(0, 15, 0, 1, 0);
        let footer = PbitFooter {
            ref_index_offset: 36,
            delta_data_offset: 36,
            sample_index_offset: 36,
            paf_data_offset: 36,
        };
        let mut buf = Vec::new();
        header.write_to(&mut buf)?;
        footer.write_to(&mut buf)?;
        std::fs::write(&path, &buf)?;

        let res = Decompressor::open(&path);
        assert!(res.is_err(), "expected rejection for zero segment_size");
        let err = res.err().unwrap().to_string();
        assert!(err.contains("segment_size"), "unexpected error: {}", err);
        Ok(())
    }

    #[test]
    fn test_decompressor_rejects_zero_kmer_len() -> Result<()> {
        // A malicious archive with kmer_len=0 must be rejected up front:
        // `open_for_append` re-segments FASTAs with this value and
        // `detect_rev_comp` calls `windows(0)`, which panics (Zero Panic).
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("bad.pbit");
        let header = PbitHeader::new(4096, 0, 18, 1, 0);
        let footer = PbitFooter {
            ref_index_offset: 36,
            delta_data_offset: 36,
            sample_index_offset: 36,
            paf_data_offset: 36,
        };
        let mut buf = Vec::new();
        header.write_to(&mut buf)?;
        footer.write_to(&mut buf)?;
        std::fs::write(&path, &buf)?;

        let res = Decompressor::open(&path);
        assert!(res.is_err(), "expected rejection for zero kmer_len");
        let err = res.err().unwrap().to_string();
        assert!(err.contains("kmer_len"), "unexpected error: {}", err);
        Ok(())
    }

    #[test]
    fn test_decompressor_rejects_huge_min_match_len() -> Result<()> {
        // A malicious archive with a huge min_match_len (and an equally huge
        // segment_size so the `min_match_len <= segment_size` check passes)
        // must be rejected up front: min_match_len drives
        // `LzDiff::prepare`'s `reference.resize(len + key_len)` padding, so an
        // unbounded value would trigger a multi-GB allocation per decode.
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("bad.pbit");
        let header = PbitHeader::new(
            MAX_PACKED_SIZE as u32 + 1,
            MAX_PACKED_SIZE as u32 + 1,
            MAX_PACKED_SIZE as u32 + 1,
            1,
            0,
        );
        let footer = PbitFooter {
            ref_index_offset: 36,
            delta_data_offset: 36,
            sample_index_offset: 36,
            paf_data_offset: 36,
        };
        let mut buf = Vec::new();
        header.write_to(&mut buf)?;
        footer.write_to(&mut buf)?;
        std::fs::write(&path, &buf)?;

        let res = Decompressor::open(&path);
        assert!(res.is_err(), "expected rejection for huge min_match_len");
        let err = res.err().unwrap().to_string();
        assert!(err.contains("min_match_len"), "unexpected error: {}", err);
        Ok(())
    }

    #[test]
    fn test_decompressor_rejects_huge_packed_size() -> Result<()> {
        // A malicious archive with an inflated delta packed_size must be
        // rejected up front (valid magic/version), not allocate a multi-GB
        // buffer in decode_delta.
        use crate::libs::fmt::twobit::write_2bit_record;
        use crate::libs::pbit::format::{write_ref_index, DeltaEncoding, DeltaMeta, RefGroupEntry};

        let dir = tempfile::tempdir()?;
        let path = dir.path().join("bad.pbit");
        let header = PbitHeader::new(4096, 15, 18, 1, 1);
        let mut buf = Vec::new();
        header.write_to(&mut buf)?; // 36 bytes

        // One reference record for "chr1" at offset 36.
        let ref_offset = buf.len() as u64;
        write_2bit_record(&mut buf, "ACGT", false)?;

        // Reference index.
        let ref_index_offset = buf.len() as u64;
        write_ref_index(
            &mut buf,
            &[RefGroupEntry {
                contig_name: "chr1".to_string(),
                ref_id: 0,
                segment_offset: ref_offset,
            }],
        )?;

        // Delta data: 1 group, 1 delta with an inflated packed_size.
        let delta_data_offset = buf.len() as u64;
        buf.extend_from_slice(&1u32.to_le_bytes()); // ref_group_count
        buf.extend_from_slice(&1u32.to_le_bytes()); // delta_count
        DeltaMeta {
            is_rev_comp: false,
            raw_length: 4,
            packed_size: MAX_PACKED_SIZE as u32 + 1,
            encoding: DeltaEncoding::LzDiff,
        }
        .write_header(&mut buf)?;

        let footer = PbitFooter {
            ref_index_offset,
            delta_data_offset,
            sample_index_offset: buf.len() as u64,
            paf_data_offset: buf.len() as u64,
        };
        footer.write_to(&mut buf)?;
        std::fs::write(&path, &buf)?;

        let res = Decompressor::open(&path);
        assert!(res.is_err(), "expected rejection for huge packed_size");
        let err = res.err().unwrap().to_string();
        assert!(err.contains("packed_size"), "unexpected error: {}", err);
        Ok(())
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
    fn test_get_sample_roundtrip_soft_mask() -> Result<()> {
        // v1005: sample soft-mask (lowercase) intervals must survive the
        // create → get_sample roundtrip losslessly (inherits 2bit semantics).
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        let ref_seq = random_dna(2000, 42);
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let sample_path = dir.path().join("sample.fa");
        let mut sample_seq = ref_seq.clone();
        // A masked (lowercase) run identical to the reference, plus some
        // uppercase mismatches elsewhere.
        for b in &mut sample_seq[500..520] {
            *b = b.to_ascii_lowercase();
        }
        sample_seq[100] = b'G';
        sample_seq[200] = b'C';
        write_fasta(sample_path.to_str().unwrap(), &[("chr1", &sample_seq)]);

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample("sample1", sample_path.to_str().unwrap())?;
        comp.finish()?;

        let mut dec = Decompressor::open(&out_path)?;
        let mut out_buf = Vec::new();
        dec.get_sample("sample1", &mut out_buf)?;
        let out_str = String::from_utf8(out_buf)?;
        let lines: Vec<&str> = out_str.lines().collect();
        assert!(lines[0].starts_with(">chr1"));
        let seq: String = lines[1..].concat();
        // Lossless: the masked run stays lowercase, mismatches stay uppercase.
        assert_eq!(seq.as_bytes(), sample_seq.as_slice());
        Ok(())
    }

    #[test]
    fn test_get_contig_roundtrip_soft_mask() -> Result<()> {
        // v1005: soft-mask (lowercase) must survive `some`/`range` extraction
        // (get_contig) losslessly, for both a full-contig and a sliced /
        // negative-strand request, matching `get_sample` (to-fa).
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        let ref_seq = random_dna(2000, 42);
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let sample_path = dir.path().join("sample.fa");
        let mut sample_seq = ref_seq.clone();
        for b in &mut sample_seq[500..520] {
            *b = b.to_ascii_lowercase();
        }
        sample_seq[100] = b'G';
        write_fasta(sample_path.to_str().unwrap(), &[("chr1", &sample_seq)]);

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample("s1", sample_path.to_str().unwrap())?;
        comp.finish()?;

        let mut dec = Decompressor::open(&out_path)?;

        // Full contig, positive strand: mask restored exactly.
        let mut full = Vec::new();
        dec.get_contig("chr1", None, None, "+", &mut full)?;
        let full_seq: String = String::from_utf8(full)?
            .lines()
            .filter(|l| !l.starts_with('>'))
            .collect();
        assert_eq!(full_seq.as_bytes(), sample_seq.as_slice());

        // Slice [490, 520): mask interval offset by the slice start.
        let mut sliced = Vec::new();
        dec.get_contig("chr1", Some(490), Some(520), "+", &mut sliced)?;
        let slice_seq: String = String::from_utf8(sliced)?
            .lines()
            .filter(|l| !l.starts_with('>'))
            .collect();
        assert_eq!(slice_seq.as_bytes(), &sample_seq[490..520]);

        // Negative strand: rev_comp preserves lowercase at transformed coords.
        let mut neg = Vec::new();
        dec.get_contig("chr1", Some(490), Some(520), "-", &mut neg)?;
        let neg_seq: String = String::from_utf8(neg)?
            .lines()
            .filter(|l| !l.starts_with('>'))
            .collect();
        let expected_neg: Vec<u8> = nt::rev_comp(&sample_seq[490..520]).collect();
        assert_eq!(neg_seq.as_bytes(), expected_neg);
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
        file.seek(SeekFrom::End(-32))?;
        let mut footer_buf = [0u8; 32];
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
        file.seek(SeekFrom::End(-32))?;
        let mut footer_buf = [0u8; 32];
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
