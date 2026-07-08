//! Compressor: builds a `.pbit` archive from reference + sample FASTA files.
//!
//! Holds a `W: Write + Seek` writer directly (no archive wrapper). The
//! reference layer is stored as standard 2bit records (reusing
//! `twobit::write_2bit_record`); sample segments are LZ-diff encoded against
//! the matching reference segment, flate2-compressed, and stored as delta
//! entries.

use anyhow::{Context, Result};
use indexmap::IndexMap;
use std::io::{Seek, SeekFrom, Write};
use std::path::Path;

use crate::libs::fmt::twobit::{read_2bit_record, write_2bit_record};
use crate::libs::nt;

use super::cigar_delta::pack_cigar;
use super::collection::Collection;
use super::decompressor::Decompressor;
use super::format::{
    read_u32_le, write_ref_index, write_u32_le, DeltaEncoding, DeltaEntry, PbitFooter, PbitHeader,
    RefGroupEntry,
};
use super::paf_index::PafQueryIndex;
use super::segment::Segment;
use crate::libs::paf::cigar::{gap_compressed_identity, CigarOp};

/// Read a FASTA file into a vector of (contig_name, sequence_bytes) pairs.
fn read_fasta(path: &str) -> Result<Vec<(String, Vec<u8>)>> {
    let mut reader = crate::libs::fmt::fa::reader(path)?;
    let mut out = Vec::new();
    for result in reader.records() {
        let record = result?;
        let name = String::from_utf8(record.name().into())?;
        let seq: Vec<u8> = record.sequence().as_ref().to_vec();
        out.push((name, seq));
    }
    Ok(out)
}

/// Split a sequence into segments of `segment_size` (last segment may be
/// shorter). Empty contigs produce no segments.
fn segment_sequence(seq: &[u8], segment_size: usize) -> Vec<&[u8]> {
    if seq.is_empty() {
        return Vec::new();
    }
    seq.chunks(segment_size).collect()
}

/// Detect orientation by sampling k-mers from `sample_seg` and checking
/// forward vs rev-comp presence in `ref_seg`. Returns `true` if rev-comp
/// gives more k-mer hits (i.e. the sample appears to be reverse-complemented
/// relative to the reference).
fn detect_rev_comp(sample_seg: &[u8], ref_seg: &[u8], kmer_len: usize) -> bool {
    if sample_seg.len() < kmer_len || ref_seg.len() < kmer_len {
        return false;
    }
    // Build a small set of k-mers sampled from the sample segment.
    let step = (sample_seg.len() / 16).max(1);
    let mut sample_kmers: Vec<&[u8]> = Vec::new();
    let mut i = 0;
    while i + kmer_len <= sample_seg.len() {
        sample_kmers.push(&sample_seg[i..i + kmer_len]);
        i += step;
    }
    // Count forward hits: how many sample k-mers appear in ref_seg forward.
    let fwd_hits = sample_kmers
        .iter()
        .filter(|k| ref_seg.windows(kmer_len).any(|w| w == **k))
        .count();
    // Count rev-comp hits: how many sample k-mers appear in rev-comp(ref_seg).
    let rc: Vec<u8> = nt::rev_comp(ref_seg).collect();
    let rc_hits = sample_kmers
        .iter()
        .filter(|k| rc.windows(kmer_len).any(|w| w == **k))
        .count();
    // Pick rev-comp if it has strictly more hits (ties go to forward).
    rc_hits > fwd_hits
}

/// Reverse-complement a byte slice into a new Vec.
fn rev_comp_vec(seq: &[u8]) -> Vec<u8> {
    nt::rev_comp(seq).collect()
}

/// Slice CIGAR to the query sub-interval [q_start, q_end) and project to the
/// target axis. Returns (sliced_ops, target_start, target_end) where
/// target_start/target_end are absolute target coordinates corresponding to
/// q_start/q_end. D ops at boundaries are excluded; D ops strictly inside
/// the segment are preserved.
fn slice_cigar_by_query(
    cigar: &[CigarOp],
    rec_qs: i32,
    rec_ts: i32,
    q_start: i32,
    q_end: i32,
) -> (Vec<CigarOp>, i32, i32) {
    let mut out: Vec<CigarOp> = Vec::new();
    let mut cur_q = rec_qs;
    let mut cur_t = rec_ts;
    let mut t_start: Option<i32> = None;
    let mut t_end: i32 = rec_ts;

    for &op in cigar {
        if cur_q >= q_end {
            break;
        }
        let qd = op.query_delta() as i32;
        let td = op.target_delta() as i32;
        let op_qs = cur_q;
        let op_qe = cur_q + qd;
        let op_ts = cur_t;

        if op.op() == 'D' {
            // D has no query span; include only if strictly inside [q_start, q_end)
            if t_start.is_some() && op_qs > q_start && op_qs < q_end {
                out.push(op);
                t_end = op_ts + td; // D advances target
            }
            cur_t = op_ts + td;
            // cur_q unchanged (qd == 0)
            continue;
        }

        // =/X/M/I: has query span
        let o_qs = op_qs.max(q_start);
        let o_qe = op_qe.min(q_end);
        if o_qe > o_qs {
            if t_start.is_none() {
                t_start = Some(match op.op() {
                    'I' => op_ts,
                    _ => op_ts + (o_qs - op_qs),
                });
            }
            t_end = match op.op() {
                'I' => op_ts,
                _ => op_ts + (o_qe - op_qs),
            };
            let overlap_len = (o_qe - o_qs) as u32;
            out.push(CigarOp::new(overlap_len, op.op()));
        }
        cur_q = op_qe;
        cur_t = op_ts + td;
    }

    let ts = t_start.unwrap_or(rec_ts);
    (out, ts, t_end)
}

/// Push or merge a CIGAR op: if the last op in `ops` has the same op char,
/// extend its length; otherwise push a new op.
fn push_or_merge(ops: &mut Vec<CigarOp>, len: u32, op_char: char) {
    match ops.last_mut() {
        Some(last) if last.op() == op_char => {
            *last = CigarOp::new(last.len() + len, op_char);
        }
        _ => ops.push(CigarOp::new(len, op_char)),
    }
}

/// Split M ops into =/X by comparing ref and sample bases, and collect X/I
/// bases into a stream. Returns (new_cigar_with_eqx, xi_bases).
fn split_m_to_eqx(
    ref_seq: &[u8],
    sample_seq: &[u8],
    cigar: &[CigarOp],
) -> Result<(Vec<CigarOp>, Vec<u8>)> {
    let mut out_ops: Vec<CigarOp> = Vec::new();
    let mut xi_bases: Vec<u8> = Vec::new();
    let mut rt: usize = 0;
    let mut si: usize = 0;

    for &op in cigar {
        let len = op.len() as usize;
        match op.op() {
            '=' => {
                push_or_merge(&mut out_ops, len as u32, '=');
                rt += len;
                si += len;
            }
            'X' => {
                if si + len > sample_seq.len() {
                    anyhow::bail!("CIGAR X exceeds sample length");
                }
                xi_bases.extend_from_slice(&sample_seq[si..si + len]);
                push_or_merge(&mut out_ops, len as u32, 'X');
                rt += len;
                si += len;
            }
            'I' => {
                if si + len > sample_seq.len() {
                    anyhow::bail!("CIGAR I exceeds sample length");
                }
                xi_bases.extend_from_slice(&sample_seq[si..si + len]);
                push_or_merge(&mut out_ops, len as u32, 'I');
                si += len;
            }
            'D' => {
                push_or_merge(&mut out_ops, len as u32, 'D');
                rt += len;
            }
            'M' => {
                if rt + len > ref_seq.len() || si + len > sample_seq.len() {
                    anyhow::bail!("CIGAR M exceeds ref/sample length");
                }
                for i in 0..len {
                    let rb = ref_seq[rt + i];
                    let sb = sample_seq[si + i];
                    if rb.eq_ignore_ascii_case(&sb) {
                        push_or_merge(&mut out_ops, 1, '=');
                    } else {
                        push_or_merge(&mut out_ops, 1, 'X');
                        xi_bases.push(sb);
                    }
                }
                rt += len;
                si += len;
            }
            other => anyhow::bail!("invalid CIGAR op: '{}'", other),
        }
    }
    Ok((out_ops, xi_bases))
}

/// Compressor: writes a `.pbit` archive.
pub struct Compressor<W: Write + Seek> {
    writer: W,
    header: PbitHeader,
    ref_groups: Vec<RefGroupEntry>,
    /// deltas[ref_group_id][delta_id] — unique deltas per ref group.
    deltas: Vec<Vec<DeltaEntry>>,
    collection: Collection,
    /// One Segment per ref_group, prepared with the (forward) reference DNA.
    segments: Vec<Segment>,
    /// Map: contig_name → Vec<ref_group_id> (reference segment indices).
    contig_ref_groups: IndexMap<String, Vec<u32>>,
    segment_size: usize,
    kmer_len: usize,
}

impl Compressor<std::io::BufWriter<std::fs::File>> {
    /// Create a new `.pbit` archive from a reference FASTA.
    ///
    /// Writes the header (placeholder offsets) + reference records (one 2bit
    /// record per segment). The caller then calls `append_sample` for each
    /// input FASTA, followed by `finish`.
    pub fn create<P: AsRef<Path>>(
        out_path: P,
        ref_fasta: &str,
        segment_size: usize,
        kmer_len: usize,
        min_match_len: u32,
    ) -> Result<Self> {
        let file = std::fs::File::create(&out_path).with_context(|| {
            format!(
                "failed to create output file: {}",
                out_path.as_ref().display()
            )
        })?;
        let writer = std::io::BufWriter::new(file);

        // Read reference FASTA and build ref_groups + segments.
        let ref_contigs = read_fasta(ref_fasta)
            .with_context(|| format!("failed to read reference FASTA: {}", ref_fasta))?;

        // We'll write the header first with a placeholder, then reference records.
        // The header's ref_records_offset is always 36 (right after the 36-byte header).
        let ref_group_count = ref_contigs
            .iter()
            .map(|(_, seq)| segment_sequence(seq, segment_size).len())
            .sum();

        let header = PbitHeader::new(
            segment_size as u32,
            kmer_len as u32,
            min_match_len,
            ref_group_count as u32,
            0, // sample_count, patched in finish()
        );

        let mut comp = Self {
            writer,
            header,
            ref_groups: Vec::new(),
            deltas: vec![Vec::new(); ref_group_count],
            collection: Collection::new(),
            segments: Vec::new(),
            contig_ref_groups: IndexMap::new(),
            segment_size,
            kmer_len,
        };

        // Write header (placeholder — ref_records_offset is already 36).
        comp.header.write_to(&mut comp.writer)?;

        // Write reference records and build the ref_groups index.
        let mut ref_group_id: u32 = 0;
        for (contig_name, seq) in &ref_contigs {
            let segs = segment_sequence(seq, segment_size);
            comp.contig_ref_groups
                .entry(contig_name.clone())
                .or_default();
            for seg in segs {
                let offset = comp.writer.stream_position()?;
                // do_mask=true preserves soft-mask (lowercase) info in 2bit record.
                let seg_str = std::str::from_utf8(seg)
                    .with_context(|| "reference segment is not valid UTF-8")?;
                write_2bit_record(&mut comp.writer, seg_str, true)?;

                let group_id = ref_group_id;
                comp.ref_groups.push(RefGroupEntry {
                    contig_name: contig_name.clone(),
                    segment_offset: offset,
                });
                comp.contig_ref_groups
                    .get_mut(contig_name)
                    .unwrap()
                    .push(group_id);

                // Prepare a Segment for this reference group.
                let mut lz = Segment::new(min_match_len);
                lz.prepare(seg);
                lz.prepare_index();
                comp.segments.push(lz);

                ref_group_id += 1;
            }
        }

        // Verify ref_group_count matches.
        debug_assert_eq!(comp.ref_groups.len() as u32, comp.header.ref_group_count);

        Ok(comp)
    }

    /// Open an existing `.pbit` for appending samples (powers `pgr pbit append`).
    /// Reads the existing header, reference records, delta data (with packed
    /// data), and collection; rebuilds Segment objects. The writer is
    /// positioned at the old ref_index_offset and the file is truncated
    /// there, ready for `append_sample` + `finish`.
    pub fn open_for_append<P: AsRef<Path>>(in_path: P) -> Result<Self> {
        let path = in_path.as_ref();

        // 1. Read archive metadata via Decompressor (opens file read-only).
        let dec = Decompressor::open(path)?;
        let header = dec.header().clone();
        let ref_groups = dec.ref_groups().to_vec();
        let collection = dec.collection_clone();
        let footer = dec.footer().clone();
        let min_match_len = header.min_match_len;
        drop(dec); // release the read-only file handle

        // 2. Reopen file for read + write.
        let file = std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .open(path)
            .with_context(|| format!("failed to open pbit file for append: {}", path.display()))?;

        // 3. Read full delta entries (header + packed_data) from delta data section.
        let mut reader = std::io::BufReader::new(file.try_clone()?);
        reader.seek(SeekFrom::Start(footer.delta_data_offset))?;
        let ref_group_count = read_u32_le(&mut reader)? as usize;
        let mut deltas: Vec<Vec<DeltaEntry>> = Vec::with_capacity(ref_group_count);
        for _ in 0..ref_group_count {
            let delta_count = read_u32_le(&mut reader)? as usize;
            let mut group = Vec::with_capacity(delta_count);
            for _ in 0..delta_count {
                group.push(DeltaEntry::read_from(&mut reader)?);
            }
            deltas.push(group);
        }

        // 4. Read reference segments and build Segment objects.
        let mut segments: Vec<Segment> = Vec::with_capacity(ref_group_count);
        let mut contig_ref_groups: IndexMap<String, Vec<u32>> = IndexMap::new();
        for (i, entry) in ref_groups.iter().enumerate() {
            reader.seek(SeekFrom::Start(entry.segment_offset))?;
            let seq = read_2bit_record(&mut reader, false, None, None, true)?;
            let seq_bytes = seq.into_bytes();
            contig_ref_groups
                .entry(entry.contig_name.clone())
                .or_default()
                .push(i as u32);
            let mut seg = Segment::new(min_match_len);
            seg.prepare(&seq_bytes);
            seg.prepare_index();
            segments.push(seg);
        }

        // 5. Truncate file at old ref_index_offset and position writer there.
        file.set_len(footer.ref_index_offset)?;
        let mut writer = std::io::BufWriter::new(file);
        writer.seek(SeekFrom::Start(footer.ref_index_offset))?;

        let segment_size = header.segment_size as usize;
        let kmer_len = header.kmer_len as usize;

        Ok(Self {
            writer,
            header,
            ref_groups,
            deltas,
            collection,
            segments,
            contig_ref_groups,
            segment_size,
            kmer_len,
        })
    }
}

impl<W: Write + Seek> Compressor<W> {
    /// Append a sample from a FASTA file. The sample name is provided by the
    /// caller (derived from the FASTA basename in the CLI layer).
    pub fn append_sample(&mut self, sample_name: &str, fasta_path: &str) -> Result<()> {
        // Ensure the sample is registered even if all contigs are unknown.
        self.collection.ensure_sample(sample_name);

        let contigs = read_fasta(fasta_path)
            .with_context(|| format!("failed to read sample FASTA: {}", fasta_path))?;

        for (contig_name, seq) in &contigs {
            // Clone to release the immutable borrow on self before calling
            // &mut self methods (encode_segment_lzdiff).
            let ref_group_ids: Vec<u32> = match self.contig_ref_groups.get(contig_name) {
                Some(ids) => ids.clone(),
                None => {
                    log::warn!(
                        "contig '{}' in sample '{}' not found in reference; skipping",
                        contig_name,
                        sample_name
                    );
                    continue;
                }
            };

            let segs = segment_sequence(seq, self.segment_size);
            if segs.is_empty() {
                // Empty contig: register with no segments.
                self.collection
                    .register_sample_contig(sample_name, contig_name);
                continue;
            }

            // Detect orientation using the first segment vs first reference segment.
            let first_ref_group = ref_group_ids[0];
            let first_ref_dna = self.segments[first_ref_group as usize].reference_dna();
            let contig_is_rev_comp = detect_rev_comp(segs[0], &first_ref_dna, self.kmer_len);

            for (seg_idx, seg) in segs.iter().enumerate() {
                self.encode_segment_lzdiff(
                    sample_name,
                    contig_name,
                    seg_idx,
                    seg,
                    &ref_group_ids,
                    contig_is_rev_comp,
                )?;
            }
        }

        Ok(())
    }

    /// LZ-diff encode one segment and append to the collection. Used by both
    /// `append_sample` and `append_sample_with_paf` (fallback path). LZ-diff
    /// segments always get `ref_start=0, ref_end=0`.
    fn encode_segment_lzdiff(
        &mut self,
        sample_name: &str,
        contig_name: &str,
        seg_idx: usize,
        seg: &[u8],
        ref_group_ids: &[u32],
        contig_is_rev_comp: bool,
    ) -> Result<()> {
        // Match to reference segment by position (clamped to last).
        let ref_idx = seg_idx.min(ref_group_ids.len() - 1);
        let ref_group_id = ref_group_ids[ref_idx];

        // Try contig-level orientation first.
        let fwd_seq: Vec<u8> = if contig_is_rev_comp {
            rev_comp_vec(seg)
        } else {
            seg.to_vec()
        };
        let fwd_delta = self.segments[ref_group_id as usize].add(&fwd_seq);
        let fwd_raw_len = fwd_seq.len() as u32;

        // If delta is large (poor match), try opposite orientation and pick smaller.
        let (delta, is_rev_comp, raw_length) = if fwd_delta.len() as u32 > fwd_raw_len / 2 {
            let alt_seq: Vec<u8> = if contig_is_rev_comp {
                seg.to_vec()
            } else {
                rev_comp_vec(seg)
            };
            let alt_delta = self.segments[ref_group_id as usize].add(&alt_seq);
            let alt_raw_len = alt_seq.len() as u32;
            if alt_delta.len() < fwd_delta.len() {
                (alt_delta, !contig_is_rev_comp, alt_raw_len)
            } else {
                (fwd_delta, contig_is_rev_comp, fwd_raw_len)
            }
        } else {
            (fwd_delta, contig_is_rev_comp, fwd_raw_len)
        };

        // flate2 compress the delta.
        let packed_data = flate2_compress(&delta)?;

        // Delta dedup: check if an identical packed_data already exists
        // in this ref_group. If so, reuse its delta_id.
        let existing = self.deltas[ref_group_id as usize]
            .iter()
            .position(|d| d.packed_data == packed_data);
        let delta_id = match existing {
            Some(id) => id as u32,
            None => {
                let entry = DeltaEntry {
                    is_rev_comp,
                    raw_length,
                    packed_data,
                    encoding: DeltaEncoding::LzDiff,
                };
                self.deltas[ref_group_id as usize].push(entry);
                (self.deltas[ref_group_id as usize].len() - 1) as u32
            }
        };

        self.collection
            .add_segment(sample_name, contig_name, ref_group_id, delta_id, 0, 0);
        Ok(())
    }

    /// Try to CIGAR-encode one segment using PAF alignments. Falls back to
    /// LZ-diff (returns Ok(false)) if: no alignment covers the segment, the
    /// best alignment doesn't fully cover it, or the target projection crosses
    /// a reference segment boundary. Returns Ok(true) if CIGAR-encoded.
    fn try_encode_segment_cigar(
        &mut self,
        sample_name: &str,
        contig_name: &str,
        seg_idx: usize,
        seg: &[u8],
        paf_index: &PafQueryIndex,
    ) -> Result<bool> {
        let seg_start = (seg_idx * self.segment_size) as i32;
        let seg_end = seg_start + seg.len() as i32;

        // 1. Look up query_id for this contig.
        let query_id = match paf_index.query_id(contig_name) {
            Some(id) => id,
            None => return Ok(false),
        };

        // 2. Query alignments overlapping [seg_start, seg_end).
        let hits = paf_index.query(query_id, seg_start, seg_end);
        if hits.is_empty() {
            return Ok(false);
        }

        // 3. Select best alignment: max coverage of [seg_start, seg_end), then max identity.
        let best = hits
            .iter()
            .max_by(|a, b| {
                let cov_a = (a.query_end.min(seg_end) - a.query_start.max(seg_start)).max(0);
                let cov_b = (b.query_end.min(seg_end) - b.query_start.max(seg_start)).max(0);
                cov_a.cmp(&cov_b).then_with(|| {
                    let id_a = gap_compressed_identity(&a.cigar);
                    let id_b = gap_compressed_identity(&b.cigar);
                    id_a.partial_cmp(&id_b).unwrap_or(std::cmp::Ordering::Equal)
                })
            })
            .unwrap();

        // 4. Check full coverage (decision 3a).
        if best.query_start > seg_start || best.query_end < seg_end {
            return Ok(false);
        }

        // 5. Slice CIGAR to [seg_start, seg_end) and project to target axis.
        let (sliced_ops, target_start, target_end) = slice_cigar_by_query(
            &best.cigar,
            best.query_start,
            best.target_start,
            seg_start,
            seg_end,
        );
        if sliced_ops.is_empty() {
            return Ok(false);
        }

        // 6. Check target doesn't cross ref segment boundary (decision 3c).
        let seg_size = self.segment_size as i32;
        let t_seg_idx_start = target_start / seg_size;
        let t_seg_idx_end = (target_end - 1) / seg_size;
        if t_seg_idx_start != t_seg_idx_end {
            return Ok(false);
        }

        // 7. Map target contig → ref_group_id.
        let ref_group_ids = match self.contig_ref_groups.get(&best.target_name) {
            Some(ids) => ids,
            None => return Ok(false),
        };
        let t_seg_idx = t_seg_idx_start as usize;
        if t_seg_idx >= ref_group_ids.len() {
            return Ok(false);
        }
        let ref_group_id = ref_group_ids[t_seg_idx];

        // 8. Compute ref_start/ref_end (relative to ref segment start).
        let ref_start = (target_start - t_seg_idx_start * seg_size) as u32;
        let ref_end = (target_end - t_seg_idx_start * seg_size) as u32;

        // 9. Get reference slice.
        let ref_dna = self.segments[ref_group_id as usize].reference_dna();
        if (ref_end as usize) > ref_dna.len() {
            return Ok(false);
        }
        let ref_slice = &ref_dna[ref_start as usize..ref_end as usize];

        // 10. Get sample slice (RC if minus strand — CIGAR describes RC(query) vs forward(target)).
        let sample_slice: Vec<u8> = if best.strand == '-' {
            rev_comp_vec(seg)
        } else {
            seg.to_vec()
        };

        // 11. Split M ops into =/X, collect X/I bases.
        let (cigar_eqx, xi_bases) = split_m_to_eqx(ref_slice, &sample_slice, &sliced_ops)?;

        // 12. Pack and store.
        let packed_data = pack_cigar(&cigar_eqx, &xi_bases)?;
        let raw_length = seg.len() as u32;
        let is_rev_comp = best.strand == '-';

        // 13. Delta dedup by packed_data.
        let existing = self.deltas[ref_group_id as usize]
            .iter()
            .position(|d| d.packed_data == packed_data);
        let delta_id = match existing {
            Some(id) => id as u32,
            None => {
                let entry = DeltaEntry {
                    is_rev_comp,
                    raw_length,
                    packed_data,
                    encoding: DeltaEncoding::Cigar,
                };
                self.deltas[ref_group_id as usize].push(entry);
                (self.deltas[ref_group_id as usize].len() - 1) as u32
            }
        };

        self.collection.add_segment(
            sample_name,
            contig_name,
            ref_group_id,
            delta_id,
            ref_start,
            ref_end,
        );
        Ok(true)
    }

    /// Append a sample using PAF-driven CIGAR encoding. Segments covered by PAF
    /// alignments are CIGAR-encoded; uncovered segments fall back to LZ-diff.
    pub fn append_sample_with_paf(
        &mut self,
        sample_name: &str,
        fasta_path: &str,
        paf_path: &str,
    ) -> Result<()> {
        self.collection.ensure_sample(sample_name);

        // Build PAF query-side index.
        let paf_index = PafQueryIndex::build_from_path(paf_path)
            .with_context(|| format!("failed to build PAF index: {}", paf_path))?;

        let contigs = read_fasta(fasta_path)
            .with_context(|| format!("failed to read sample FASTA: {}", fasta_path))?;

        for (contig_name, seq) in &contigs {
            // Clone to release the immutable borrow on self before calling
            // &mut self methods (try_encode_segment_cigar / encode_segment_lzdiff).
            let ref_group_ids: Vec<u32> = match self.contig_ref_groups.get(contig_name) {
                Some(ids) => ids.clone(),
                None => {
                    log::warn!(
                        "contig '{}' in sample '{}' not found in reference; skipping",
                        contig_name,
                        sample_name
                    );
                    continue;
                }
            };

            let segs = segment_sequence(seq, self.segment_size);
            if segs.is_empty() {
                self.collection
                    .register_sample_contig(sample_name, contig_name);
                continue;
            }

            // Detect orientation (hint for LZ-diff fallback only).
            let first_ref_group = ref_group_ids[0];
            let first_ref_dna = self.segments[first_ref_group as usize].reference_dna();
            let contig_is_rev_comp = detect_rev_comp(segs[0], &first_ref_dna, self.kmer_len);

            for (seg_idx, seg) in segs.iter().enumerate() {
                // Try CIGAR encoding first; fall back to LZ-diff.
                let encoded = self.try_encode_segment_cigar(
                    sample_name,
                    contig_name,
                    seg_idx,
                    seg,
                    &paf_index,
                )?;
                if !encoded {
                    self.encode_segment_lzdiff(
                        sample_name,
                        contig_name,
                        seg_idx,
                        seg,
                        &ref_group_ids,
                        contig_is_rev_comp,
                    )?;
                }
            }
        }
        Ok(())
    }

    /// Finalize: write Reference Index → Delta Data → Sample Index → Footer →
    /// patch Header sample_count. Consumes the compressor.
    pub fn finish(mut self) -> Result<()> {
        // Patch header sample_count.
        self.header.sample_count = self.collection.sample_count() as u32;

        // Seek to the end of reference records (current writer position).
        let ref_index_offset = self.writer.stream_position()?;

        // Write Reference Index.
        write_ref_index(&mut self.writer, &self.ref_groups)?;

        // Write Delta Data.
        let delta_data_offset = self.writer.stream_position()?;
        write_u32_le(&mut self.writer, self.deltas.len() as u32)?;
        for group_deltas in &self.deltas {
            write_u32_le(&mut self.writer, group_deltas.len() as u32)?;
            for entry in group_deltas {
                entry.write_to(&mut self.writer)?;
            }
        }

        // Write Sample Index (collection, flate2-compressed).
        let sample_index_offset = self.writer.stream_position()?;
        let collection_bytes = self.collection.serialize()?;
        self.writer.write_all(&collection_bytes)?;

        // Write Footer.
        let footer = PbitFooter {
            ref_index_offset,
            delta_data_offset,
            sample_index_offset,
        };
        footer.write_to(&mut self.writer)?;

        // Patch header (sample_count may have changed; rewrite at offset 0).
        self.writer.seek(SeekFrom::Start(0))?;
        self.header.write_to(&mut self.writer)?;

        self.writer.flush()?;
        Ok(())
    }

    /// Set the command line string stored in the collection.
    pub fn set_cmd_line(&mut self, cmd: &str) {
        self.collection.cmd_line = cmd.to_string();
    }
}

/// flate2-compress a byte slice.
pub(crate) fn flate2_compress(data: &[u8]) -> Result<Vec<u8>> {
    use std::io::Write;
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(data)?;
    Ok(encoder.finish()?)
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_create_and_finish_empty() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", b"ACGTACGTACGTACGT")]);
        let out_path = dir.path().join("out.pbit");

        let comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.finish()?;

        assert!(out_path.exists());
        Ok(())
    }

    #[test]
    fn test_create_with_one_sample() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        let ref_seq = random_dna(5000, 42);
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let sample_path = dir.path().join("sample.fa");
        let mut sample_seq = ref_seq.clone();
        // Introduce a few SNPs.
        sample_seq[100] = b'G';
        sample_seq[200] = b'C';
        sample_seq[300] = b'T';
        write_fasta(sample_path.to_str().unwrap(), &[("chr1", &sample_seq)]);

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample("sample1", sample_path.to_str().unwrap())?;
        comp.finish()?;

        // Verify the file is non-empty and starts with the magic.
        let mut file = std::fs::File::open(&out_path)?;
        let header = PbitHeader::read_from(&mut file)?;
        assert_eq!(header.magic, super::super::format::PBIT_MAGIC);
        assert_eq!(header.sample_count, 1);
        assert_eq!(header.ref_group_count, 2); // 5000 bp / 4096 = 2 segments
        Ok(())
    }

    #[test]
    fn test_create_multiple_samples_dedup() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        let ref_seq = random_dna(2000, 42);
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        // Two identical samples → delta dedup should collapse them.
        let s1_path = dir.path().join("s1.fa");
        let s2_path = dir.path().join("s2.fa");
        write_fasta(s1_path.to_str().unwrap(), &[("chr1", &ref_seq)]);
        write_fasta(s2_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample("s1", s1_path.to_str().unwrap())?;
        comp.append_sample("s2", s2_path.to_str().unwrap())?;
        comp.finish()?;

        // Read back and verify.
        let mut file = std::fs::File::open(&out_path)?;
        let header = PbitHeader::read_from(&mut file)?;
        assert_eq!(header.sample_count, 2);
        assert_eq!(header.ref_group_count, 1); // 2000 bp < 4096 → 1 segment
        Ok(())
    }

    #[test]
    fn test_skip_unknown_contig() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        let ref_seq = random_dna(1000, 42);
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let sample_path = dir.path().join("sample.fa");
        let sample_seq = random_dna(1000, 99);
        write_fasta(
            sample_path.to_str().unwrap(),
            &[("unknown_contig", &sample_seq)],
        );

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample("sample1", sample_path.to_str().unwrap())?;
        comp.finish()?;

        // The sample should have 0 contigs (all skipped).
        let mut file = std::fs::File::open(&out_path)?;
        let header = PbitHeader::read_from(&mut file)?;
        assert_eq!(header.sample_count, 1);
        Ok(())
    }

    #[test]
    fn test_detect_rev_comp() {
        // Non-palindromic reference (not equal to its own rev-comp).
        let ref_seq = b"AAATCGGGCTAGCCATAGGCCGATTAAGCCGA";
        let sample_fwd = ref_seq;
        let sample_rev: Vec<u8> = nt::rev_comp(ref_seq).collect();
        // Forward sample should not trigger rev-comp.
        assert!(!detect_rev_comp(sample_fwd, ref_seq, 8));
        // Rev-comp sample should trigger rev-comp.
        assert!(detect_rev_comp(&sample_rev, ref_seq, 8));
    }

    #[test]
    fn test_segment_sequence() {
        let seq = vec![b'A'; 10];
        let segs = segment_sequence(&seq, 4);
        assert_eq!(segs.len(), 3);
        assert_eq!(segs[0].len(), 4);
        assert_eq!(segs[1].len(), 4);
        assert_eq!(segs[2].len(), 2);

        // Empty sequence → no segments.
        assert!(segment_sequence(&[], 4).is_empty());
    }

    #[test]
    fn test_flate2_roundtrip() -> Result<()> {
        let data = b"hello world hello world hello world";
        let compressed = flate2_compress(data)?;
        use std::io::Read;
        let mut decoder = flate2::read::GzDecoder::new(&compressed[..]);
        let mut decompressed = Vec::new();
        decoder.read_to_end(&mut decompressed)?;
        assert_eq!(decompressed, data);
        Ok(())
    }

    #[test]
    fn test_open_for_append() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        let ref_seq = random_dna(2000, 42);
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let s1_path = dir.path().join("s1.fa");
        let s1_seq = introduce_snps(&ref_seq, 100);
        write_fasta(s1_path.to_str().unwrap(), &[("chr1", &s1_seq)]);

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample("s1", s1_path.to_str().unwrap())?;
        comp.finish()?;

        // Append a second sample.
        let s2_path = dir.path().join("s2.fa");
        let s2_seq = introduce_snps(&ref_seq, 200);
        write_fasta(s2_path.to_str().unwrap(), &[("chr1", &s2_seq)]);

        let mut comp = Compressor::open_for_append(&out_path)?;
        comp.append_sample("s2", s2_path.to_str().unwrap())?;
        comp.finish()?;

        // Verify both samples are present and extract correctly.
        let mut dec = crate::libs::pbit::decompressor::Decompressor::open(&out_path)?;
        assert_eq!(dec.list_samples(), vec!["s1", "s2"]);

        let mut buf = Vec::new();
        dec.get_sample("s2", &mut buf)?;
        let out_str = String::from_utf8(buf)?;
        let lines: Vec<&str> = out_str.lines().collect();
        let seq: String = lines[1..].concat();
        let expected =
            String::from_utf8(s2_seq.iter().map(|&c| c.to_ascii_uppercase()).collect()).unwrap();
        assert_eq!(seq, expected);
        Ok(())
    }

    /// Introduce SNPs at every 100th position (helper for append test).
    fn introduce_snps(seq: &[u8], seed: u64) -> Vec<u8> {
        use rand::rngs::StdRng;
        use rand::Rng;
        use rand::SeedableRng;
        let mut rng = StdRng::seed_from_u64(seed);
        let mut out = seq.to_vec();
        for i in (0..out.len()).step_by(100) {
            out[i] = match out[i] {
                b'A' => {
                    if rng.random_range(0u8..3) == 0 {
                        b'C'
                    } else {
                        b'G'
                    }
                }
                _ => b'A',
            };
        }
        out
    }

    /// Build a single PAF line string.
    #[allow(clippy::too_many_arguments)]
    fn paf_line(
        qname: &str,
        qlen: u32,
        qs: u32,
        qe: u32,
        strand: &str,
        tname: &str,
        tlen: u32,
        ts: u32,
        te: u32,
        cigar: &str,
    ) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t0\t255\tcg:Z:{}",
            qname, qlen, qs, qe, strand, tname, tlen, ts, te, cigar
        )
    }

    /// Write a PAF file with one line per string in `lines`.
    fn write_paf(path: &str, lines: &[String]) {
        use std::io::Write;
        let mut f = std::fs::File::create(path).unwrap();
        for line in lines {
            writeln!(f, "{}", line).unwrap();
        }
    }

    /// Extract the sequence (concatenated, sans header) from a FASTA buffer.
    fn extract_fasta_seq(buf: &[u8]) -> String {
        let s = String::from_utf8_lossy(buf);
        s.lines().skip(1).collect::<String>().trim().to_string()
    }

    // ── slice_cigar_by_query tests ───────────────────────────

    #[test]
    fn test_slice_cigar_by_query_pure_match() {
        // CIGAR: 100=, rec_qs=0, rec_ts=0. Slice [20, 50).
        let ops = crate::libs::paf::cigar::parse_cigar("100=").unwrap();
        let (sliced, ts, te) = slice_cigar_by_query(&ops, 0, 0, 20, 50);
        assert_eq!(sliced.len(), 1);
        assert_eq!(sliced[0].op(), '=');
        assert_eq!(sliced[0].len(), 30);
        assert_eq!(ts, 20);
        assert_eq!(te, 50);
    }

    #[test]
    fn test_slice_cigar_by_query_with_indel() {
        // CIGAR: 10=5I10=5D10=, rec_qs=0, rec_ts=0. Slice [5, 30).
        // Trace:
        //   10=: q[0,10)  t[0,10)  → overlap q[5,10) t[5,10)  → 5=,  t_start=5,  t_end=10
        //   5I: q[10,15) t[10,10) → overlap q[10,15) t[10,10) → 5I,  t_end=10
        //   10=: q[15,25) t[10,20)→ overlap q[15,25) t[10,20) → 10=, t_end=20
        //   5D: q[25,25) t[20,25) → D inside (25>5 && 25<30)  → 5D,  t_end=25
        //   10=: q[25,35) t[25,35)→ overlap q[25,30) t[25,30) → 5=,  t_end=30
        let ops = crate::libs::paf::cigar::parse_cigar("10=5I10=5D10=").unwrap();
        let (sliced, ts, te) = slice_cigar_by_query(&ops, 0, 0, 5, 30);
        assert_eq!(sliced.len(), 5);
        assert_eq!(sliced[0], CigarOp::new(5, '='));
        assert_eq!(sliced[1], CigarOp::new(5, 'I'));
        assert_eq!(sliced[2], CigarOp::new(10, '='));
        assert_eq!(sliced[3], CigarOp::new(5, 'D'));
        assert_eq!(sliced[4], CigarOp::new(5, '='));
        assert_eq!(ts, 5);
        assert_eq!(te, 30);
    }

    // ── split_m_to_eqx tests ─────────────────────────────────

    #[test]
    fn test_split_m_to_eqx_all_match() {
        // ref == sample, pure M CIGAR → all become =.
        let ref_seq = b"ACGTACGT";
        let sample_seq = b"ACGTACGT";
        let cigar = crate::libs::paf::cigar::parse_cigar("8M").unwrap();
        let (ops, xi) = split_m_to_eqx(ref_seq, sample_seq, &cigar).unwrap();
        assert_eq!(ops.len(), 1);
        assert_eq!(ops[0], CigarOp::new(8, '='));
        assert!(xi.is_empty());
    }

    #[test]
    fn test_split_m_to_eqx_with_mismatches() {
        // ref = ACGTACGT, sample = ACGAACGT, cigar = 8M
        // Position 3: ref=T, sample=A → X, xi=[A]
        let ref_seq = b"ACGTACGT";
        let sample_seq = b"ACGAACGT";
        let cigar = crate::libs::paf::cigar::parse_cigar("8M").unwrap();
        let (ops, xi) = split_m_to_eqx(ref_seq, sample_seq, &cigar).unwrap();
        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0], CigarOp::new(3, '='));
        assert_eq!(ops[1], CigarOp::new(1, 'X'));
        assert_eq!(ops[2], CigarOp::new(4, '='));
        assert_eq!(xi, vec![b'A']);
    }

    // ── append_sample_with_paf tests ─────────────────────────

    #[test]
    fn test_append_sample_with_paf_plus_strand() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_seq = random_dna(2000, 42);
        let ref_path = dir.path().join("ref.fa");
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        // Sample = ref with 3 SNPs at 100, 200, 300.
        let mut sample_seq = ref_seq.clone();
        sample_seq[100] = b'G';
        sample_seq[200] = b'C';
        sample_seq[300] = b'T';
        let sample_path = dir.path().join("sample.fa");
        write_fasta(sample_path.to_str().unwrap(), &[("chr1", &sample_seq)]);

        // PAF: + strand, full coverage, CIGAR describes the SNPs.
        let paf_path = dir.path().join("sample.paf");
        write_paf(
            paf_path.to_str().unwrap(),
            &[paf_line(
                "chr1",
                2000,
                0,
                2000,
                "+",
                "chr1",
                2000,
                0,
                2000,
                "100=1X99=1X99=1X1699=",
            )],
        );

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample_with_paf(
            "sample1",
            sample_path.to_str().unwrap(),
            paf_path.to_str().unwrap(),
        )?;
        comp.finish()?;

        // Decompress and verify.
        let mut dec = crate::libs::pbit::decompressor::Decompressor::open(&out_path)?;
        let mut buf = Vec::new();
        dec.get_sample("sample1", &mut buf)?;
        let got = extract_fasta_seq(&buf);
        let expected = String::from_utf8(sample_seq.clone()).unwrap();
        assert_eq!(got, expected);
        Ok(())
    }

    #[test]
    fn test_append_sample_with_paf_minus_strand() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_seq = random_dna(2000, 42);
        let ref_path = dir.path().join("ref.fa");
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        // Sample = RC(ref). PAF says - strand, so CIGAR describes RC(sample) vs ref = ref vs ref.
        let sample_seq: Vec<u8> = nt::rev_comp(&ref_seq).collect();
        let sample_path = dir.path().join("sample.fa");
        write_fasta(sample_path.to_str().unwrap(), &[("chr1", &sample_seq)]);

        let paf_path = dir.path().join("sample.paf");
        write_paf(
            paf_path.to_str().unwrap(),
            &[paf_line(
                "chr1", 2000, 0, 2000, "-", "chr1", 2000, 0, 2000, "2000=",
            )],
        );

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample_with_paf(
            "sample1",
            sample_path.to_str().unwrap(),
            paf_path.to_str().unwrap(),
        )?;
        comp.finish()?;

        let mut dec = crate::libs::pbit::decompressor::Decompressor::open(&out_path)?;
        let mut buf = Vec::new();
        dec.get_sample("sample1", &mut buf)?;
        let got = extract_fasta_seq(&buf);
        let expected = String::from_utf8(sample_seq.clone()).unwrap();
        assert_eq!(got, expected);
        Ok(())
    }

    #[test]
    fn test_append_sample_with_paf_partial_coverage() -> Result<()> {
        // ref = 5000 bp → 2 segments (4096 + 904). PAF only covers first segment.
        let dir = tempfile::tempdir()?;
        let ref_seq = random_dna(5000, 42);
        let ref_path = dir.path().join("ref.fa");
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let sample_seq = ref_seq.clone();
        let sample_path = dir.path().join("sample.fa");
        write_fasta(sample_path.to_str().unwrap(), &[("chr1", &sample_seq)]);

        // PAF covers only [0, 4096) — second segment falls back to LZ-diff.
        let paf_path = dir.path().join("sample.paf");
        write_paf(
            paf_path.to_str().unwrap(),
            &[paf_line(
                "chr1", 5000, 0, 4096, "+", "chr1", 5000, 0, 4096, "4096=",
            )],
        );

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample_with_paf(
            "sample1",
            sample_path.to_str().unwrap(),
            paf_path.to_str().unwrap(),
        )?;
        comp.finish()?;

        let mut dec = crate::libs::pbit::decompressor::Decompressor::open(&out_path)?;
        let mut buf = Vec::new();
        dec.get_sample("sample1", &mut buf)?;
        let got = extract_fasta_seq(&buf);
        let expected = String::from_utf8(sample_seq.clone()).unwrap();
        assert_eq!(got, expected);
        Ok(())
    }

    #[test]
    fn test_append_sample_with_paf_empty_paf() -> Result<()> {
        // Empty PAF file → all segments fall back to LZ-diff.
        let dir = tempfile::tempdir()?;
        let ref_seq = random_dna(2000, 42);
        let ref_path = dir.path().join("ref.fa");
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let mut sample_seq = ref_seq.clone();
        sample_seq[100] = b'G';
        let sample_path = dir.path().join("sample.fa");
        write_fasta(sample_path.to_str().unwrap(), &[("chr1", &sample_seq)]);

        let paf_path = dir.path().join("empty.paf");
        write_paf(paf_path.to_str().unwrap(), &[]);

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample_with_paf(
            "sample1",
            sample_path.to_str().unwrap(),
            paf_path.to_str().unwrap(),
        )?;
        comp.finish()?;

        let mut dec = crate::libs::pbit::decompressor::Decompressor::open(&out_path)?;
        let mut buf = Vec::new();
        dec.get_sample("sample1", &mut buf)?;
        let got = extract_fasta_seq(&buf);
        let expected = String::from_utf8(sample_seq.clone()).unwrap();
        assert_eq!(got, expected);
        Ok(())
    }
}
