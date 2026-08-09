//! Compressor: builds a `.pbit` archive from reference + sample FASTA files.
//!
//! Holds a `W: Write + Seek` writer directly (no archive wrapper). The
//! reference layer is stored as standard 2bit records (reusing
//! `twobit::write_2bit_record`); sample segments are LZ-diff encoded against
//! the matching reference segment, flate2-compressed, and stored as delta
//! entries.

use anyhow::{bail, Context, Result};
use indexmap::IndexMap;
use std::collections::HashMap;
use std::io::{Seek, SeekFrom, Write};
use std::path::Path;

use crate::libs::fmt::twobit::{read_2bit_record, write_2bit_record};
use crate::libs::nt;

use super::cigar_delta::pack_cigar;
use super::collection::Collection;
use super::decompressor::Decompressor;
use super::format::{
    read_u32_le, write_ref_index, write_ref_table, write_string, write_u32_le, DeltaEncoding,
    DeltaEntry, PafRecovery, PbitFooter, PbitHeader, RefGroupEntry, RefTableEntry,
};
use super::paf_index::PafQueryIndex;
use super::segment::Segment;
use crate::libs::paf::cigar::{extract_cigar, gap_compressed_identity, CigarOp};

/// Read a FASTA file into a vector of (contig_name, sequence_bytes) pairs.
fn read_fasta(path: &str) -> Result<Vec<(String, Vec<u8>)>> {
    let mut reader = crate::libs::fmt::seq::SeqReader::new(path)?;
    let mut rec = crate::libs::fmt::seq::SeqRecord::new();
    let mut out = Vec::new();
    while reader.read_record(&mut rec)? {
        let name = String::from_utf8(rec.name().to_vec())?;
        let seq: Vec<u8> = rec.sequence().to_vec();
        out.push((name, seq));
    }
    Ok(out)
}

/// Extract soft-mask intervals (lowercase runs) as `(start, size)` in 0-based
/// coordinates, same semantics as 2bit `mask_blocks` (v1005).
fn extract_mask_blocks(seq: &[u8]) -> Vec<(u32, u32)> {
    let mut blocks = Vec::new();
    let mut i = 0usize;
    while i < seq.len() {
        if seq[i].is_ascii_lowercase() {
            let start = i;
            while i < seq.len() && seq[i].is_ascii_lowercase() {
                i += 1;
            }
            blocks.push((start as u32, (i - start) as u32));
        } else {
            i += 1;
        }
    }
    blocks
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
///
/// Orientation is detected once per contig using the first segment. This
/// assumption holds for typical whole-chromosome alignments but may be
/// incorrect for contigs with internal inversions or rearrangements; such
/// segments fall back to LZ-diff's opposite-orientation retry logic.
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

/// Convert a forward-strand query sub-interval [seg_start, seg_end) into
/// CIGAR coordinate space for a '-' strand record. CIGAR for '-' records
/// describes RC(query) vs target, so CIGAR position 0 corresponds to
/// forward position query_end - 1. Returns `(rc_start, rc_end)`.
fn forward_to_rc_coords(seg_start: i32, seg_end: i32, query_end: i32) -> (i32, i32) {
    (query_end - seg_end, query_end - seg_start)
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
                if rt + len > ref_seq.len() || si + len > sample_seq.len() {
                    bail!("CIGAR '=' exceeds ref/sample length");
                }
                push_or_merge(&mut out_ops, len as u32, '=');
                rt += len;
                si += len;
            }
            'X' => {
                if rt + len > ref_seq.len() {
                    bail!("CIGAR X exceeds ref length");
                }
                if si + len > sample_seq.len() {
                    bail!("CIGAR X exceeds sample length");
                }
                xi_bases.extend_from_slice(&sample_seq[si..si + len]);
                push_or_merge(&mut out_ops, len as u32, 'X');
                rt += len;
                si += len;
            }
            'I' => {
                if si + len > sample_seq.len() {
                    bail!("CIGAR I exceeds sample length");
                }
                xi_bases.extend_from_slice(&sample_seq[si..si + len]);
                push_or_merge(&mut out_ops, len as u32, 'I');
                si += len;
            }
            'D' => {
                if rt + len > ref_seq.len() {
                    bail!("CIGAR D exceeds ref length");
                }
                push_or_merge(&mut out_ops, len as u32, 'D');
                rt += len;
            }
            'M' => {
                if rt + len > ref_seq.len() || si + len > sample_seq.len() {
                    bail!("CIGAR M exceeds ref/sample length");
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
            other => bail!("invalid CIGAR op: '{}'", other),
        }
    }
    if rt != ref_seq.len() || si != sample_seq.len() {
        bail!(
            "CIGAR consumed ref={}/{} sample={}/{}",
            rt,
            ref_seq.len(),
            si,
            sample_seq.len()
        );
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
    /// ref_group_id → global reference offset of the segment's first base
    /// (v1007 CIGAR deltas reference reference-file-global coordinates).
    ref_seg_starts: Vec<u64>,
    collection: Collection,
    /// One Segment per ref_group, prepared with the (forward) reference DNA.
    segments: Vec<Segment>,
    /// Map: contig_name → Vec<ref_group_id> (reference segment indices).
    contig_ref_groups: IndexMap<String, Vec<u32>>,
    /// Lazy canonical-k-mer → reference segment index for content-based
    /// matching when sample contig names do not match reference contigs.
    ref_kmer_index: Option<HashMap<u64, Vec<u32>>>,
    /// Per-reference metadata (name, group range, embedded-index offsets).
    ref_meta: Vec<RefTableEntry>,
    /// Reference a sample routes to during `append_sample` (set per sample).
    cur_ref_id: u32,
    /// PAF recovery data per sample (v1009): (sample name, big-chain ms
    /// table (record_id, ms), verbatim small-chain PAF rows).
    paf_data: Vec<PafRecovery>,
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
        Self::create_multi(
            out_path,
            &[ref_fasta],
            segment_size,
            kmer_len,
            min_match_len,
        )
    }

    /// Create a new `.pbit` archive from one or more reference FASTA files.
    ///
    /// Each reference genome gets a distinct `ref_id` and its own segment
    /// group range; samples route to one reference (see `append_sample`).
    pub fn create_multi<P: AsRef<Path>>(
        out_path: P,
        ref_fastas: &[&str],
        segment_size: usize,
        kmer_len: usize,
        min_match_len: u32,
    ) -> Result<Self> {
        anyhow::ensure!(
            !ref_fastas.is_empty(),
            "at least one reference FASTA is required"
        );
        let file = std::fs::File::create(&out_path).with_context(|| {
            format!(
                "failed to create output file: {}",
                out_path.as_ref().display()
            )
        })?;
        let writer = std::io::BufWriter::new(file);

        // We'll write the header first with a placeholder, then reference records.
        // The header's ref_records_offset is always 36 (right after the 36-byte header).
        let mut ref_group_count = 0usize;
        let mut all_ref_contigs: Vec<Vec<(String, Vec<u8>)>> = Vec::with_capacity(ref_fastas.len());
        for ref_fasta in ref_fastas {
            let ref_contigs = read_fasta(ref_fasta)
                .with_context(|| format!("failed to read reference FASTA: {}", ref_fasta))?;
            ref_group_count += ref_contigs
                .iter()
                .map(|(_, seq)| segment_sequence(seq, segment_size).len())
                .sum::<usize>();
            all_ref_contigs.push(ref_contigs);
        }

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
            ref_seg_starts: Vec::with_capacity(ref_group_count),
            collection: Collection::new(),
            segments: Vec::new(),
            contig_ref_groups: IndexMap::new(),
            ref_kmer_index: None,
            ref_meta: Vec::new(),
            cur_ref_id: 0,
            paf_data: Vec::new(),
            segment_size,
            kmer_len,
        };

        // Write header (placeholder — ref_records_offset is already 36).
        comp.header.write_to(&mut comp.writer)?;

        // Write reference records and build the ref_groups index.
        let mut ref_group_id: u32 = 0;
        let mut ref_pos = 0u64;
        for (ref_id, ref_contigs) in all_ref_contigs.iter().enumerate() {
            let group_start = ref_group_id;
            let mut group_count = 0u32;
            for (contig_name, seq) in ref_contigs {
                let segs = segment_sequence(seq, segment_size);
                let groups = comp
                    .contig_ref_groups
                    .entry(contig_name.clone())
                    .or_default();
                for seg in segs {
                    let offset = comp.writer.stream_position()?;
                    // do_mask=true preserves soft-mask (lowercase) info in 2bit record.
                    let seg_str = std::str::from_utf8(seg)
                        .with_context(|| "reference segment is not valid UTF-8")?;
                    write_2bit_record(&mut comp.writer, seg_str, true)?;
                    comp.ref_seg_starts.push(ref_pos);
                    ref_pos += seg.len() as u64;

                    let group_id = ref_group_id;
                    comp.ref_groups.push(RefGroupEntry {
                        contig_name: contig_name.clone(),
                        ref_id: ref_id as u32,
                        segment_offset: offset,
                    });
                    groups.push(group_id);

                    // Prepare a Segment for this reference group.
                    let mut lz = Segment::new(min_match_len);
                    lz.prepare(seg);
                    lz.prepare_index();
                    comp.segments.push(lz);

                    ref_group_id += 1;
                    group_count += 1;
                }
            }
            let ref_name = crate::libs::io::get_basename(ref_fastas[ref_id])
                .unwrap_or_else(|| ref_fastas[ref_id].to_string());
            comp.ref_meta.push(RefTableEntry {
                ref_name,
                group_start,
                group_count,
            });
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
        anyhow::ensure!(
            header.segment_size <= i32::MAX as u32,
            "archive segment_size {} exceeds i32::MAX; archive is corrupt or malicious",
            header.segment_size
        );
        let ref_groups = dec.ref_groups().to_vec();
        let ref_seg_starts = dec.ref_seg_starts().to_vec();
        let collection = dec.collection_clone();
        let footer = dec.footer().clone();
        let dec_ref_table = dec.ref_table().to_vec();
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
            ref_seg_starts,
            collection,
            segments,
            contig_ref_groups,
            ref_kmer_index: None,
            ref_meta: dec_ref_table.clone(),
            cur_ref_id: 0,
            paf_data: Vec::new(),
            segment_size,
            kmer_len,
        })
    }
}

impl<W: Write + Seek> Compressor<W> {
    /// Set the reference (ref_id) that the next `append_sample` call routes
    /// to. Defaults to 0.
    pub fn set_cur_ref_id(&mut self, ref_id: u32) {
        self.cur_ref_id = ref_id;
    }

    /// Names of the reference genomes in the archive, in ref_id order.
    pub fn ref_names(&self) -> Vec<&str> {
        self.ref_meta.iter().map(|r| r.ref_name.as_str()).collect()
    }

    /// Canonical 2-bit k-mer (min of k-mer and its reverse complement) as u64.
    fn canonical_kmer(seq: &[u8], pos: usize, k: usize) -> Option<u64> {
        if pos + k > seq.len() {
            return None;
        }
        let mut mer = 0u64;
        for &b in &seq[pos..pos + k] {
            let v = match b.to_ascii_uppercase() {
                b'A' => 0u64,
                b'C' => 1u64,
                b'G' => 2u64,
                b'T' => 3u64,
                _ => return None,
            };
            mer = (mer << 2) | v;
        }
        let mut rc = 0u64;
        let mut m = mer;
        for _ in 0..k {
            rc = (rc << 2) | (3 - (m & 3));
            m >>= 2;
        }
        Some(mer.min(rc))
    }

    /// Build (once) the reference canonical-k-mer → segment index used by the
    /// content-based LZ-diff fallback.
    fn ensure_ref_kmer_index(&mut self) {
        if self.ref_kmer_index.is_some() {
            return;
        }
        let k = self.kmer_len;
        let mut idx: HashMap<u64, Vec<u32>> = HashMap::new();
        for (gid, seg) in self.segments.iter().enumerate() {
            let dna = seg.reference_dna();
            for pos in 0..dna.len().saturating_sub(k - 1) {
                if let Some(km) = Self::canonical_kmer(&dna, pos, k) {
                    idx.entry(km).or_default().push(gid as u32);
                }
            }
        }
        self.ref_kmer_index = Some(idx);
    }

    /// Reference segment with the most shared canonical k-mers with `seg`.
    fn best_ref_group(&mut self, seg: &[u8]) -> Option<u32> {
        self.ensure_ref_kmer_index();
        let idx = self.ref_kmer_index.as_ref()?;
        let k = self.kmer_len;
        let mut votes: HashMap<u32, u32> = HashMap::new();
        for pos in 0..seg.len().saturating_sub(k - 1) {
            if let Some(km) = Self::canonical_kmer(seg, pos, k) {
                if let Some(gids) = idx.get(&km) {
                    for &g in gids {
                        *votes.entry(g).or_default() += 1;
                    }
                }
            }
        }
        votes.into_iter().max_by_key(|&(_, v)| v).map(|(g, _)| g)
    }

    /// Whether a sample `name` already exists in the archive. Appending a
    /// sample whose name collides with an existing one would silently merge
    /// their segments (corrupting the sample on extract); callers should
    /// reject this before calling `append_sample`.
    pub fn has_sample(&self, name: &str) -> bool {
        self.collection.samples.contains_key(name)
    }

    /// Append a sample from a FASTA file. The sample name is provided by the
    /// caller (derived from the FASTA basename in the CLI layer).
    pub fn append_sample(&mut self, sample_name: &str, fasta_path: &str) -> Result<()> {
        // Ensure the sample is registered even if all contigs are unknown.
        self.collection.ensure_sample(sample_name);

        let contigs = read_fasta(fasta_path)
            .with_context(|| format!("failed to read sample FASTA: {}", fasta_path))?;

        for (contig_name, seq) in &contigs {
            // Soft-mask intervals come from the original case-preserving
            // sequence; encoding uses the uppercase copy (2bit semantics).
            let mask_blocks = extract_mask_blocks(seq);
            let seq_upper: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
            let segs = segment_sequence(&seq_upper, self.segment_size);
            if segs.is_empty() {
                // Empty contig: register with no segments.
                self.collection
                    .register_sample_contig(sample_name, contig_name);
                continue;
            }
            let Some(meta) = self.ref_meta.get(self.cur_ref_id as usize) else {
                anyhow::bail!(
                    "invalid reference id {} ({} references)",
                    self.cur_ref_id,
                    self.ref_meta.len()
                );
            };
            let (group_start, group_count) = (meta.group_start, meta.group_count);
            let ref_group_ids: Vec<u32> = self
                .contig_ref_groups
                .get(contig_name)
                .map(|ids| {
                    ids.iter()
                        .copied()
                        .filter(|id| *id >= group_start && *id < group_start + group_count)
                        .collect()
                })
                .unwrap_or_default();
            if ref_group_ids.is_empty() {
                // Content-based fallback (design §8.5 route 1): sample contig
                // names differ from reference names; match each segment to its
                // best-matching reference segment by canonical k-mer overlap
                // and LZ-encode against it (no contig-name identity required).
                let mut any = false;
                for (seg_idx, seg) in segs.iter().enumerate() {
                    if let Some(gid) = self.best_ref_group(seg) {
                        self.encode_segment_lzdiff(
                            sample_name,
                            contig_name,
                            seg_idx,
                            (seg_idx as u32) * (self.segment_size as u32),
                            seg,
                            &[gid],
                            false,
                        )?;
                    } else {
                        // No content match in the reference: store the segment
                        // verbatim so the archive stays lossless.
                        self.encode_segment_raw(
                            sample_name,
                            contig_name,
                            (seg_idx as u32) * (self.segment_size as u32),
                            seg,
                        )?;
                    }
                    any = true;
                }
                if any {
                    self.collection
                        .register_sample_contig(sample_name, contig_name)
                        .mask_blocks = mask_blocks;
                }
                continue;
            }
            // Register the contig entry with its soft-mask intervals (v1005);
            // subsequent `add_segment` calls reuse this entry.
            {
                let cs = self
                    .collection
                    .register_sample_contig(sample_name, contig_name);
                cs.mask_blocks = mask_blocks;
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
                    (seg_idx as u32) * (self.segment_size as u32),
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
    #[allow(clippy::too_many_arguments)]
    fn encode_segment_lzdiff(
        &mut self,
        sample_name: &str,
        contig_name: &str,
        seg_idx: usize,
        q_start: u32,
        seg: &[u8],
        ref_group_ids: &[u32],
        contig_is_rev_comp: bool,
    ) -> Result<()> {
        anyhow::ensure!(
            !ref_group_ids.is_empty(),
            "encode_segment_lzdiff: no reference groups for contig '{}' in sample '{}'",
            contig_name,
            sample_name
        );
        // Match to reference segment by position (clamped to last). For a
        // reverse-complemented contig, the sample's segment i is the rev-comp
        // of reference segment N-1-i (the reference's last segments appear
        // first in the RC sample), so route in reverse order; routing forward
        // would still decode correctly but yield a poor LZ-diff match.
        let last = ref_group_ids.len().saturating_sub(1);
        let ref_idx = if contig_is_rev_comp {
            last.saturating_sub(seg_idx)
        } else {
            seg_idx.min(last)
        };
        let ref_group_id = ref_group_ids[ref_idx];

        // Try contig-level orientation first.
        let fwd_seq: Vec<u8> = if contig_is_rev_comp {
            rev_comp_vec(seg)
        } else {
            seg.to_vec()
        };
        let fwd_delta = self.segments[ref_group_id as usize].add(&fwd_seq)?;
        let fwd_raw_len = fwd_seq.len() as u32;

        // If delta is large (poor match), try opposite orientation and pick smaller.
        let (delta, is_rev_comp, raw_length) = if fwd_delta.len() as u32 > fwd_raw_len / 2 {
            let alt_seq: Vec<u8> = if contig_is_rev_comp {
                seg.to_vec()
            } else {
                rev_comp_vec(seg)
            };
            let alt_delta = self.segments[ref_group_id as usize].add(&alt_seq)?;
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

        // Delta dedup: check if an identical packed_data with the same
        // orientation already exists in this ref_group. is_rev_comp is part
        // of the delta header and cannot be shared across orientations
        // (e.g., an empty delta for forward vs rev-comp reference).
        let existing = self.deltas[ref_group_id as usize]
            .iter()
            .position(|d| d.packed_data == packed_data && d.is_rev_comp == is_rev_comp);
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

        self.collection.add_segment(
            sample_name,
            contig_name,
            ref_group_id,
            delta_id,
            0,
            0,
            q_start,
            u32::MAX,
        );
        Ok(())
    }

    /// Store a segment verbatim as a standard 2bit record (v1008+) so
    /// sequences with no matching reference content are preserved
    /// losslessly. 2bit is ~7 pp smaller than flate2(ASCII) for DNA and
    /// needs no secondary compression. Raw deltas are attached to ref_group
    /// 0 and their reference segment is never read during decode.
    fn encode_segment_raw(
        &mut self,
        sample_name: &str,
        contig_name: &str,
        q_start: u32,
        seg: &[u8],
    ) -> Result<()> {
        let seg_str = std::str::from_utf8(seg).with_context(|| "raw segment is not valid UTF-8")?;
        let mut buf = Vec::with_capacity(seg.len() / 2);
        {
            let mut cursor = std::io::Cursor::new(&mut buf);
            write_2bit_record(&mut cursor, seg_str, false)?;
        }
        let packed_data = buf;
        let gid = 0u32;
        let existing = self.deltas[gid as usize]
            .iter()
            .position(|d| d.packed_data == packed_data && d.encoding == DeltaEncoding::Raw);
        let delta_id = match existing {
            Some(id) => id as u32,
            None => {
                let entry = DeltaEntry {
                    is_rev_comp: false,
                    raw_length: seg.len() as u32,
                    packed_data,
                    encoding: DeltaEncoding::Raw,
                };
                self.deltas[gid as usize].push(entry);
                (self.deltas[gid as usize].len() - 1) as u32
            }
        };
        self.collection.add_segment(
            sample_name,
            contig_name,
            gid,
            delta_id,
            0,
            0,
            q_start,
            u32::MAX,
        );
        Ok(())
    }

    /// Read an arbitrary reference interval in reference-file-global
    /// coordinates from the in-memory reference segments (mirror of
    /// `Decompressor::read_ref_interval`).
    fn read_ref_interval_local(&self, start: u64, end: u64) -> Result<Vec<u8>> {
        if start >= end {
            return Ok(Vec::new());
        }
        let mut idx = match self.ref_seg_starts.binary_search(&start) {
            Ok(i) => i,
            Err(0) => anyhow::bail!("read_ref_interval_local: start {start} before first segment"),
            Err(i) => i - 1,
        };
        let mut out = Vec::with_capacity((end - start) as usize);
        let mut cur = start;
        while cur < end {
            let seg_start = self.ref_seg_starts[idx];
            let seg_end = if idx + 1 < self.ref_seg_starts.len() {
                self.ref_seg_starts[idx + 1]
            } else {
                seg_start + self.segments[idx].reference_dna().len() as u64
            };
            let seg = self.segments[idx].reference_dna();
            let lo = (cur - seg_start) as usize;
            let hi = ((end - seg_start) as usize).min((seg_end - seg_start) as usize);
            out.extend_from_slice(&seg[lo..hi]);
            cur = seg_start + hi as u64;
            idx += 1;
        }
        Ok(out)
    }

    /// Select main (primary) chains for one sample: chain-level greedy by
    /// query-covered segment count (desc) → identity (desc) → input order.
    /// A chain whose query interval overlaps any already-chosen main chain
    /// is a small (secondary) chain; its PAF row is stored verbatim so
    /// `to-paf` can reproduce it (2026-08-09, replaces the v1009 length
    /// threshold). Records without `cg:Z` are never main chains.
    fn select_main_chains(paf_index: &PafQueryIndex, segment_size: usize) -> Vec<bool> {
        let seg_size = segment_size.max(1) as i64;
        let mut candidates: Vec<(u32, u32, f64, u32)> = Vec::new(); // (qs, qe, gi, record_id)
        for (id, line) in paf_index.records.iter().enumerate() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 12 {
                continue;
            }
            let tag_strs: Vec<String> = f[12..].iter().map(|s| s.to_string()).collect();
            if !tag_strs.iter().any(|t| t.starts_with("cg:Z:")) {
                continue;
            }
            let (Ok(qs), Ok(qe)) = (f[2].parse::<u32>(), f[3].parse::<u32>()) else {
                continue;
            };
            let Ok(cigar) = extract_cigar(&tag_strs) else {
                continue;
            };
            if cigar.is_empty() {
                continue;
            }
            let gi = gap_compressed_identity(&cigar);
            candidates.push((qs, qe, gi, id as u32));
        }
        // Sort: covered segment count desc → identity desc → input order.
        candidates.sort_by(|a, b| {
            let segs_a = ((a.1 - a.0) as i64 + seg_size - 1) / seg_size;
            let segs_b = ((b.1 - b.0) as i64 + seg_size - 1) / seg_size;
            segs_b
                .cmp(&segs_a)
                .then_with(|| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal))
                .then_with(|| a.3.cmp(&b.3))
        });
        let mut main = vec![false; paf_index.records.len()];
        let mut chosen: Vec<(u32, u32)> = Vec::new();
        for (qs, qe, _, id) in candidates {
            if chosen.iter().any(|&(s, e)| s < qe && qs < e) {
                continue; // overlaps a chosen main chain → small chain
            }
            main[id as usize] = true;
            chosen.push((qs, qe));
        }
        main
    }

    /// Try to CIGAR-encode the PAF-covered sub-interval of one segment.
    /// Returns `Ok(Some((q0, q1)))` with the covered sample-contig interval
    /// when a CIGAR delta was stored, or `Ok(None)` when the segment has no
    /// PAF coverage. Since v1007 the delta may reference any reference
    /// interval (global coordinates), so partial coverage and target spans
    /// across reference blocks are both supported; the caller encodes the
    /// uncovered parts with LZ-diff/Raw. Any chain with a cg:Z tag may
    /// encode a segment (2026-08-09); main/small selection only decides whether
    /// the chain is rebuilt or stored verbatim in the PAF recovery data.
    fn try_encode_segment_cigar(
        &mut self,
        sample_name: &str,
        contig_name: &str,
        seg_idx: usize,
        seg: &[u8],
        paf_index: &PafQueryIndex,
    ) -> Result<Option<(u32, u32, u32)>> {
        let seg_start = (seg_idx * self.segment_size) as i32;
        let seg_end = seg_start + seg.len() as i32;

        // 1. Look up query_id for this contig.
        let query_id = match paf_index.query_id(contig_name) {
            Some(id) => id,
            None => return Ok(None),
        };

        // 2. Query alignments overlapping [seg_start, seg_end).
        let hits = paf_index.query(query_id, seg_start, seg_end);
        if hits.is_empty() {
            return Ok(None);
        }

        // 3. Select the chain covering this segment most (2026-08-09: any chain
        // with cg:Z may encode; main/small only affects PAF recovery).
        // If no chain covers the segment, fall back to LZ-diff/Raw.
        let cover_hits: Vec<_> = hits.iter().collect();
        if cover_hits.is_empty() {
            return Ok(None);
        }
        let best = cover_hits
            .into_iter()
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
        let q0 = best.query_start.max(seg_start);
        let q1 = best.query_end.min(seg_end);
        if q0 >= q1 {
            return Ok(None);
        }

        // 5. Slice CIGAR to the segment interval and project to the target
        // axis. For '-' strand the CIGAR describes RC(query) vs forward
        // (target): CIGAR op 0 corresponds to RC(query) position 0 (= forward
        // query_end-1), and the CIGAR traverses forward query coords from
        // high to low. Convert forward [seg_start, seg_end) to RC coords
        // [query_end - seg_end, query_end - seg_start) and slice with
        // rec_qs = 0 (CIGAR origin = RC(query) position 0). For '+' strand
        // the CIGAR traverses forward query coords low→high, so forward
        // coords are used directly with rec_qs = query_start.
        let (sliced_ops, target_start, target_end) = if best.strand == '+' {
            slice_cigar_by_query(&best.cigar, best.query_start, best.target_start, q0, q1)
        } else {
            let (rc_start, rc_end) = forward_to_rc_coords(q0, q1, best.query_end);
            slice_cigar_by_query(&best.cigar, 0, best.target_start, rc_start, rc_end)
        };
        if sliced_ops.is_empty() {
            return Ok(None);
        }

        // 6. Map target contig → ref_group_id (the block containing the
        // target start; the delta may still span further blocks, which is
        // fine since ref coordinates are global since v1007).
        let Some(meta) = self.ref_meta.get(self.cur_ref_id as usize) else {
            return Ok(None);
        };
        let ref_group_ids: Vec<u32> = match self.contig_ref_groups.get(&best.target_name) {
            Some(ids) => ids
                .iter()
                .copied()
                .filter(|id| *id >= meta.group_start && *id < meta.group_start + meta.group_count)
                .collect(),
            None => return Ok(None),
        };
        let seg_size = self.segment_size as i32;
        let t_seg_idx_start = (target_start / seg_size).max(0);
        let t_seg_idx = t_seg_idx_start as usize;
        if t_seg_idx >= ref_group_ids.len() {
            return Ok(None);
        }
        let ref_group_id = ref_group_ids[t_seg_idx];

        // 7. Compute ref_start/ref_end in reference-file-global coordinates.
        let ref_start = self.ref_seg_starts[ref_group_id as usize]
            + (target_start - t_seg_idx_start * seg_size) as u64;
        let ref_end = self.ref_seg_starts[ref_group_id as usize]
            + (target_end - t_seg_idx_start * seg_size) as u64;

        // 8. Get reference slice (global coordinates → segment slices).
        let ref_slice = self.read_ref_interval_local(ref_start, ref_end)?;

        // 9. Get sample slice for the covered interval (RC if minus strand —
        // CIGAR describes RC(query) vs forward(target)).
        let seg_lo = (q0 - seg_start) as usize;
        let seg_hi = (q1 - seg_start) as usize;
        let sample_slice: Vec<u8> = if best.strand == '-' {
            let rc: Vec<u8> = rev_comp_vec(&seg[seg_lo..seg_hi]);
            rc
        } else {
            seg[seg_lo..seg_hi].to_vec()
        };

        // 11. Split M ops into =/X, collect X/I bases.
        let (cigar_eqx, xi_bases) = split_m_to_eqx(&ref_slice, &sample_slice, &sliced_ops)?;

        // 12. Pack and store.
        let raw_length = (q1 - q0) as u32;
        let is_rev_comp = best.strand == '-';
        // A segment identical to its reference interval (all '=' ops, no
        // X/I/D) is stored as a zero-payload pointer to the interval
        // (AGC-style Identity, v1010): the segment carries ref_start/ref_end,
        // and all identity segments sharing orientation + length reuse one
        // delta entry. Otherwise pack the CIGAR delta as before.
        let is_identity = cigar_eqx.iter().all(|op| op.op() == '=');
        let (packed_data, encoding) = if is_identity {
            (Vec::new(), DeltaEncoding::Identity)
        } else {
            (pack_cigar(&cigar_eqx, &xi_bases)?, DeltaEncoding::Cigar)
        };

        // 13. Delta dedup: identity by encoding + orientation + length (the
        // interval lives in the segment), CIGAR by payload + orientation.
        let existing = if is_identity {
            self.deltas[ref_group_id as usize].iter().position(|d| {
                d.encoding == DeltaEncoding::Identity
                    && d.is_rev_comp == is_rev_comp
                    && d.raw_length == raw_length
            })
        } else {
            self.deltas[ref_group_id as usize]
                .iter()
                .position(|d| d.packed_data == packed_data && d.is_rev_comp == is_rev_comp)
        };
        let delta_id = match existing {
            Some(id) => id as u32,
            None => {
                let entry = DeltaEntry {
                    is_rev_comp,
                    raw_length,
                    packed_data,
                    encoding,
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
            ref_start as u32,
            ref_end as u32,
            q0 as u32,
            best.record_id,
        );
        Ok(Some((q0 as u32, q1 as u32, best.record_id)))
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
        // Chain-level main/small selection (2026-08-09): main chains are encoded
        // and rebuilt; small chains (overlapping a chosen main chain, or
        // without cg:Z) are stored verbatim for PAF recovery.
        let main_chains = Self::select_main_chains(&paf_index, self.segment_size);

        let contigs = read_fasta(fasta_path)
            .with_context(|| format!("failed to read sample FASTA: {}", fasta_path))?;

        // Resolve the current reference's group range so the LZ-diff fallback
        // (below) routes segments to the correct reference. Without this, a
        // contig name shared across references would fall back against the
        // wrong reference's segment in multi-reference archives.
        // Copy the reference's group range into locals (u32 is Copy) so the
        // filter closure below does not borrow `self` immutably while later
        // calls borrow `self` mutably.
        let (ref_group_start, ref_group_count) = match self.ref_meta.get(self.cur_ref_id as usize) {
            Some(meta) => (meta.group_start, meta.group_count),
            None => {
                anyhow::bail!(
                    "invalid reference id {} ({} references)",
                    self.cur_ref_id,
                    self.ref_meta.len()
                );
            }
        };

        // PAF record ids used by CIGAR encoding across this sample's contigs
        // (v1009 recovery data: big chains get their ms stored, the rest are
        // stored verbatim so `to-paf` can reproduce the original PAF).
        let mut used_records: std::collections::HashSet<u32> = std::collections::HashSet::new();
        let mut used_count: std::collections::HashMap<u32, u32> = std::collections::HashMap::new();

        for (contig_name, seq) in &contigs {
            // Soft-mask intervals come from the original case-preserving
            // sequence; encoding uses the uppercase copy (2bit semantics).
            let mask_blocks = extract_mask_blocks(seq);
            let seq_upper: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
            let segs = segment_sequence(&seq_upper, self.segment_size);
            if segs.is_empty() {
                self.collection
                    .register_sample_contig(sample_name, contig_name);
                continue;
            }

            // Reference groups matched by contig name (LZ-diff fallback only;
            // the CIGAR path maps sample contigs via the PAF index and does
            // not require name identity between sample and reference).
            let ref_group_ids: Vec<u32> = self
                .contig_ref_groups
                .get(contig_name)
                .map(|ids| {
                    ids.iter()
                        .copied()
                        .filter(|id| {
                            *id >= ref_group_start && *id < ref_group_start + ref_group_count
                        })
                        .collect()
                })
                .unwrap_or_default();

            // Orientation hint for the LZ-diff fallback (only meaningful when
            // the contig name matches a reference contig).
            let contig_is_rev_comp = if ref_group_ids.is_empty() {
                false
            } else {
                let first_ref_group = ref_group_ids[0];
                let first_ref_dna = self.segments[first_ref_group as usize].reference_dna();
                detect_rev_comp(segs[0], &first_ref_dna, self.kmer_len)
            };

            let mut any_encoded = false;
            for (seg_idx, seg) in segs.iter().enumerate() {
                let seg_start = (seg_idx as u32) * (self.segment_size as u32);
                let seg_end = seg_start + seg.len() as u32;
                // Try CIGAR encoding of the PAF-covered sub-interval first.
                let covered = self.try_encode_segment_cigar(
                    sample_name,
                    contig_name,
                    seg_idx,
                    seg,
                    &paf_index,
                )?;
                if let Some((q0, q1, record_id)) = covered {
                    any_encoded = true;
                    used_records.insert(record_id);
                    *used_count.entry(record_id).or_default() += 1;
                    // Encode the uncovered parts: prefer LZ-diff content
                    // matching (AGC-style, much smaller than raw), fall back
                    // to verbatim storage so the contig stays gapless and
                    // lossless.
                    if seg_start < q0 {
                        let sub = &seg[0..(q0 - seg_start) as usize];
                        if let Some(gid) = self.best_ref_group(sub) {
                            self.encode_segment_lzdiff(
                                sample_name,
                                contig_name,
                                0,
                                seg_start,
                                sub,
                                &[gid],
                                false,
                            )?;
                        } else {
                            self.encode_segment_raw(sample_name, contig_name, seg_start, sub)?;
                        }
                    }
                    if q1 < seg_end {
                        let sub = &seg[(q1 - seg_start) as usize..];
                        if let Some(gid) = self.best_ref_group(sub) {
                            self.encode_segment_lzdiff(
                                sample_name,
                                contig_name,
                                0,
                                q1,
                                sub,
                                &[gid],
                                false,
                            )?;
                        } else {
                            self.encode_segment_raw(sample_name, contig_name, q1, sub)?;
                        }
                    }
                    continue;
                }
                // No PAF coverage for this segment: LZ-diff / Raw fallback.
                if !ref_group_ids.is_empty() {
                    self.encode_segment_lzdiff(
                        sample_name,
                        contig_name,
                        seg_idx,
                        (seg_idx as u32) * (self.segment_size as u32),
                        seg,
                        &ref_group_ids,
                        contig_is_rev_comp,
                    )?;
                    any_encoded = true;
                } else if let Some(gid) = self.best_ref_group(seg) {
                    // Content-based fallback: no contig-name match and no
                    // PAF coverage; encode against the best-matching
                    // reference segment by canonical k-mer overlap.
                    self.encode_segment_lzdiff(
                        sample_name,
                        contig_name,
                        seg_idx,
                        (seg_idx as u32) * (self.segment_size as u32),
                        seg,
                        &[gid],
                        false,
                    )?;
                    any_encoded = true;
                } else {
                    // No PAF coverage and no content match: store the segment
                    // verbatim so the archive stays lossless.
                    self.encode_segment_raw(sample_name, contig_name, seg_start, seg)?;
                    any_encoded = true;
                }
            }

            if any_encoded {
                // Attach soft-mask intervals after encoding (segments were
                // registered on demand by `add_segment`).
                self.collection
                    .register_sample_contig(sample_name, contig_name)
                    .mask_blocks = mask_blocks;
            }
        }

        // v1009: PAF recovery data for this sample — a "big chain" (span >=
        // threshold) that was CIGAR-encoded across ALL its segments is
        // rebuilt at chain level (only its ms is stored); every other record
        // (small chains — overlapping a chosen main chain or without cg:Z —
        // and main chains that lost some segment to encoding failure) is
        // stored verbatim so `to-paf` can reproduce the original PAF exactly.
        let seg_size = self.segment_size as i32;
        let mut big_ms: Vec<(u32, i32)> = Vec::new();
        let mut small: Vec<String> = Vec::new();
        for (i, line) in paf_index.records.iter().enumerate() {
            let id = i as u32;
            let f: Vec<&str> = line.split('\t').collect();
            let coords = f
                .get(2)
                .and_then(|s| s.parse::<i32>().ok())
                .zip(f.get(3).and_then(|s| s.parse::<i32>().ok()));
            let Some((qs, qe)) = coords else {
                continue;
            };
            let span_segs = if qs < qe {
                (qe - 1) / seg_size - qs / seg_size + 1
            } else {
                0
            };
            let encoded = used_count.get(&id).copied().unwrap_or(0);
            let is_main = main_chains.get(id as usize).copied().unwrap_or(false);
            let is_complete_big = is_main && encoded >= span_segs.max(1) as u32;
            if is_complete_big {
                big_ms.push((id, paf_index.record_ms(id).unwrap_or(0)));
            } else {
                small.push(line.clone());
            }
        }
        big_ms.sort_unstable();
        self.paf_data.push((sample_name.to_string(), big_ms, small));
        Ok(())
    }

    /// Finalize: write Reference Index → Delta Data → Sample Index → Footer →
    /// patch Header sample_count. Consumes the compressor.
    pub fn finish(mut self) -> Result<()> {
        // Patch header sample_count.
        self.header.sample_count = self.collection.sample_count() as u32;

        // Seek to the end of reference records (current writer position).
        let ref_index_offset = self.writer.stream_position()?;

        // Write Reference Index: per-group entries + reference table.
        write_ref_index(&mut self.writer, &self.ref_groups)?;
        write_ref_table(&mut self.writer, &self.ref_meta)?;

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

        // Write PAF recovery data (v1009): verbatim small-chain PAF rows plus
        // the ms table of CIGAR-encoded records, so `pbit to-paf` can
        // reproduce the original PAF.
        let paf_data_offset = self.writer.stream_position()?;
        self.write_paf_data()?;

        // Write Footer.
        let footer = PbitFooter {
            ref_index_offset,
            delta_data_offset,
            sample_index_offset,
            paf_data_offset,
        };
        footer.write_to(&mut self.writer)?;

        // Patch header (sample_count may have changed; rewrite at offset 0).
        self.writer.seek(SeekFrom::Start(0))?;
        self.header.write_to(&mut self.writer)?;

        self.writer.flush()?;
        Ok(())
    }

    /// Append a new reference genome at the current writer position (which
    /// must be the truncation point before the Reference Index, i.e. after
    /// `open_for_append`). Register it in `ref_meta`.
    pub fn append_reference(&mut self, ref_fasta: &str) -> Result<()> {
        let ref_id = self.ref_meta.len() as u32;
        let ref_contigs = read_fasta(ref_fasta)
            .with_context(|| format!("failed to read reference FASTA: {}", ref_fasta))?;
        let group_start = self.ref_groups.len() as u32;
        let mut group_count = 0u32;
        let mut ref_pos = self
            .ref_seg_starts
            .last()
            .copied()
            .map(|s| {
                s + self.segments[self.ref_seg_starts.len() - 1]
                    .reference_dna()
                    .len() as u64
            })
            .unwrap_or(0);
        for (contig_name, seq) in &ref_contigs {
            let segs = segment_sequence(seq, self.segment_size);
            let groups = self
                .contig_ref_groups
                .entry(contig_name.clone())
                .or_default();
            for seg in segs {
                let offset = self.writer.stream_position()?;
                let seg_str = std::str::from_utf8(seg)
                    .with_context(|| "reference segment is not valid UTF-8")?;
                write_2bit_record(&mut self.writer, seg_str, true)?;
                self.ref_seg_starts.push(ref_pos);
                ref_pos += seg.len() as u64;
                let group_id = self.ref_groups.len() as u32;
                self.ref_groups.push(RefGroupEntry {
                    contig_name: contig_name.clone(),
                    ref_id,
                    segment_offset: offset,
                });
                groups.push(group_id);
                let mut lz = Segment::new(self.header.min_match_len);
                lz.prepare(seg);
                lz.prepare_index();
                self.segments.push(lz);
                group_count += 1;
            }
        }
        self.deltas.resize(self.ref_groups.len(), Vec::new());
        self.header.ref_group_count = self.ref_groups.len() as u32;
        let ref_name =
            crate::libs::io::get_basename(ref_fasta).unwrap_or_else(|| ref_fasta.to_string());
        self.ref_meta.push(RefTableEntry {
            ref_name,
            group_start,
            group_count,
        });
        Ok(())
    }

    /// Set the command line string stored in the collection.
    pub fn set_cmd_line(&mut self, cmd: &str) {
        self.collection.cmd_line = cmd.to_string();
    }

    /// Serialize the PAF recovery data (v1009) and append it to the writer:
    /// per sample, the (record_id, ms) table of CIGAR-encoded chains and the
    /// verbatim PAF rows of all other chains, flate2-compressed.
    fn write_paf_data(&mut self) -> Result<()> {
        let mut raw = Vec::new();
        write_u32_le(&mut raw, self.paf_data.len() as u32)?;
        for (sample, big_ms, small) in &self.paf_data {
            write_string(&mut raw, sample)?;
            write_u32_le(&mut raw, big_ms.len() as u32)?;
            for (record_id, ms) in big_ms {
                write_u32_le(&mut raw, *record_id)?;
                write_u32_le(&mut raw, *ms as u32)?;
            }
            write_u32_le(&mut raw, small.len() as u32)?;
            for line in small {
                write_string(&mut raw, line)?;
            }
        }
        let compressed = flate2_compress(&raw)?;
        self.writer.write_all(&compressed)?;
        Ok(())
    }
}

/// flate2-compress a byte slice.
fn flate2_compress(data: &[u8]) -> Result<Vec<u8>> {
    Ok(crate::libs::bgzf::gzip_compress(data, 6)?)
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
    fn test_create_sets_cmd_line() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", b"ACGTACGTACGTACGT")]);
        let out_path = dir.path().join("out.pbit");

        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.set_cmd_line("pgr pbit create -r ref.fa -o out.pbit -s 4096 -k 15 -l 18");
        comp.finish()?;

        let dec = Decompressor::open(&out_path)?;
        assert!(
            dec.collection().cmd_line.contains("pgr pbit create"),
            "cmd_line should record create command: {}",
            dec.collection().cmd_line
        );
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
    fn test_raw_fallback_lossless() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        let ref_seq = random_dna(1000, 42);
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        // A contig with no content match in the reference must be stored
        // verbatim (Raw delta) and round-trip exactly, not silently skipped.
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

        let mut file = std::fs::File::open(&out_path)?;
        let header = PbitHeader::read_from(&mut file)?;
        assert_eq!(header.sample_count, 1);

        let mut dec = crate::libs::pbit::decompressor::Decompressor::open(&out_path)?;
        let mut buf = Vec::new();
        dec.get_sample("sample1", &mut buf)?;
        let out_str = String::from_utf8(buf)?;
        let lines: Vec<&str> = out_str.lines().collect();
        let seq: String = lines[1..].concat();
        assert_eq!(seq, String::from_utf8(sample_seq).unwrap());
        Ok(())
    }

    #[test]
    fn test_append_rev_comp_sample_multi_segment() -> Result<()> {
        // A reverse-complemented sample against a multi-segment reference must
        // route its segments in reverse order (segment i ↔ ref segment N-1-i).
        // Roundtrip must still reproduce the sample exactly (uppercase).
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        let ref_seq = random_dna(5000, 42); // 2 segments (4096 + 904)
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let mut sample_fwd = ref_seq.clone();
        sample_fwd[100] = b'G';
        sample_fwd[4600] = b'C';
        let sample_seq: Vec<u8> = nt::rev_comp(&sample_fwd).collect();
        let sample_path = dir.path().join("sample.fa");
        write_fasta(sample_path.to_str().unwrap(), &[("chr1", &sample_seq)]);

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample("sample1", sample_path.to_str().unwrap())?;
        comp.finish()?;

        let mut dec = crate::libs::pbit::decompressor::Decompressor::open(&out_path)?;
        let mut buf = Vec::new();
        dec.get_sample("sample1", &mut buf)?;
        let out_str = String::from_utf8(buf)?;
        let lines: Vec<&str> = out_str.lines().collect();
        let seq: String = lines[1..].concat();
        let expected =
            String::from_utf8(sample_seq.iter().map(|&c| c.to_ascii_uppercase()).collect())
                .unwrap();
        assert_eq!(seq, expected);

        // The sample spans 2 reference segments; both must round-trip.
        let header =
            crate::libs::pbit::format::PbitHeader::read_from(&mut std::fs::File::open(&out_path)?)?;
        assert_eq!(header.ref_group_count, 2);
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
        let decompressed = crate::libs::bgzf::gzip_decompress(&compressed, data.len())?;
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
        comp.set_cmd_line("pgr pbit append out.pbit -o out.pbit");
        comp.append_sample("s2", s2_path.to_str().unwrap())?;
        comp.finish()?;

        // Verify both samples are present and extract correctly.
        let mut dec = crate::libs::pbit::decompressor::Decompressor::open(&out_path)?;
        assert!(
            dec.collection().cmd_line.contains("pgr pbit append"),
            "cmd_line should record append command: {}",
            dec.collection().cmd_line
        );
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

    #[test]
    fn test_split_m_to_eqx_error_eq_overflow() {
        // `=` op longer than ref/sample must be rejected (previously the `=`
        // branch advanced cursors without bounds checks, unlike X/I/M).
        let ref_seq = b"ACGTACGT";
        let sample_seq = b"ACGTACGT";
        let cigar = crate::libs::paf::cigar::parse_cigar("10=").unwrap();
        let err = split_m_to_eqx(ref_seq, sample_seq, &cigar).unwrap_err();
        assert!(err.to_string().contains("exceeds ref/sample length"));
    }

    #[test]
    fn test_split_m_to_eqx_error_underconsumed() {
        // CIGAR that does not fully consume ref/sample must be rejected rather
        // than silently returning ops/bases that don't cover the sequences.
        let ref_seq = b"ACGTACGT";
        let sample_seq = b"ACGTACGT";
        let cigar = crate::libs::paf::cigar::parse_cigar("4=").unwrap();
        let err = split_m_to_eqx(ref_seq, sample_seq, &cigar).unwrap_err();
        assert!(err.to_string().contains("consumed ref=4/8 sample=4/8"));
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
    fn test_append_sample_with_paf_cross_assembly_names() -> Result<()> {
        // Regression: sample contigs named differently from the reference
        // (draft-assembly vs reference naming) must still be encodable via
        // the PAF/CIGAR path. The old code gated every contig on a
        // reference-contig-name lookup before trying CIGAR, so cross-assembly
        // samples were silently skipped.
        let dir = tempfile::tempdir()?;
        let ref_seq = random_dna(2000, 42);
        let ref_path = dir.path().join("ref.fa");
        write_fasta(ref_path.to_str().unwrap(), &[("refc", &ref_seq)]);

        // Sample = ref with 3 SNPs, under a different contig name.
        let mut sample_seq = ref_seq.clone();
        sample_seq[100] = b'G';
        sample_seq[200] = b'C';
        sample_seq[300] = b'T';
        let sample_path = dir.path().join("sample.fa");
        write_fasta(sample_path.to_str().unwrap(), &[("sampc", &sample_seq)]);

        let paf_path = dir.path().join("sample.paf");
        write_paf(
            paf_path.to_str().unwrap(),
            &[paf_line(
                "sampc",
                2000,
                0,
                2000,
                "+",
                "refc",
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
    fn test_append_sample_with_paf_multi_segment_gap_free() -> Result<()> {
        // 12000 bp reference = 3 segments (4096 + 4096 + 3808). A single
        // full-coverage CIGAR without internal gaps must encode all segments.
        let dir = tempfile::tempdir()?;
        let ref_seq = random_dna(12000, 42);
        let ref_path = dir.path().join("ref.fa");
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let mut sample_seq = ref_seq.clone();
        sample_seq[100] = b'G';
        sample_seq[5000] = b'C';
        sample_seq[9000] = b'T';
        let sample_path = dir.path().join("sample.fa");
        write_fasta(sample_path.to_str().unwrap(), &[("sampc", &sample_seq)]);

        let paf_path = dir.path().join("sample.paf");
        write_paf(
            paf_path.to_str().unwrap(),
            &[paf_line(
                "sampc",
                12000,
                0,
                12000,
                "+",
                "chr1",
                12000,
                0,
                12000,
                "100=1X4899=1X3999=1X3999=",
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
    fn test_append_sample_with_paf_deletion_crosses_segment() -> Result<()> {
        // A deletion whose target interval crosses a 4096-bp reference segment
        // boundary is rejected by the CIGAR path (design constraint 3), but
        // the LZ content fallback recovers the segment losslessly against the
        // best-matching reference segment, so the sample round-trips fully.
        let dir = tempfile::tempdir()?;
        let ref_seq = random_dna(12000, 42);
        let ref_path = dir.path().join("ref.fa");
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        // Sample = ref with 100 bp deleted at position 5000 (target interval
        // [5000, 5100) lies within ref segment 1, but the query segment
        // [4096, 8192) maps to target [4096, 8192+100) crossing the boundary).
        let mut sample_seq = Vec::with_capacity(11900);
        sample_seq.extend_from_slice(&ref_seq[..5000]);
        sample_seq.extend_from_slice(&ref_seq[5100..]);
        let sample_path = dir.path().join("sample.fa");
        write_fasta(sample_path.to_str().unwrap(), &[("sampc", &sample_seq)]);

        let paf_path = dir.path().join("sample.paf");
        write_paf(
            paf_path.to_str().unwrap(),
            &[paf_line(
                "sampc",
                11900,
                0,
                11900,
                "+",
                "chr1",
                12000,
                0,
                12000,
                "5000=100D6900=",
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
    fn test_append_sample_with_paf_indel_breaks_phase() -> Result<()> {
        // A 1-bp insertion early in the sample shifts the phase of every
        // downstream segment: the CIGAR path requires phase-aligned target
        // intervals and drops those segments, but the LZ content fallback
        // recovers them losslessly against the best-matching reference
        // segment, so the sample round-trips fully.
        let dir = tempfile::tempdir()?;
        let ref_seq = random_dna(12000, 42);
        let ref_path = dir.path().join("ref.fa");
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        // Sample = ref with one base inserted at position 100.
        let mut sample_seq = Vec::with_capacity(12001);
        sample_seq.extend_from_slice(&ref_seq[..100]);
        sample_seq.push(b'G');
        sample_seq.extend_from_slice(&ref_seq[100..]);
        let sample_path = dir.path().join("sample.fa");
        write_fasta(sample_path.to_str().unwrap(), &[("sampc", &sample_seq)]);

        let paf_path = dir.path().join("sample.paf");
        write_paf(
            paf_path.to_str().unwrap(),
            &[paf_line(
                "sampc",
                12001,
                0,
                12001,
                "+",
                "chr1",
                12000,
                0,
                12000,
                "100=1I11900=",
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
    fn test_append_sample_content_match_cross_assembly() -> Result<()> {
        // Sample contig names differ from the reference and no PAF is
        // provided: the LZ-diff content fallback must match each segment to
        // its best reference segment by canonical k-mer overlap and round-trip
        // losslessly (design §8.5 route 1).
        let dir = tempfile::tempdir()?;
        let ref_seq = random_dna(12000, 42);
        let ref_path = dir.path().join("ref.fa");
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        let mut sample_seq = ref_seq.clone();
        sample_seq[100] = b'G';
        sample_seq[5000] = b'C';
        sample_seq[9000] = b'T';
        let sample_path = dir.path().join("sample.fa");
        write_fasta(sample_path.to_str().unwrap(), &[("sampc", &sample_seq)]);

        let out_path = dir.path().join("out.pbit");
        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        comp.append_sample("sample1", sample_path.to_str().unwrap())?;
        comp.finish()?;

        let mut dec = crate::libs::pbit::decompressor::Decompressor::open(&out_path)?;
        let mut buf = Vec::new();
        dec.get_sample("sample1", &mut buf)?;
        let got = extract_fasta_seq(&buf);
        let expected =
            String::from_utf8(sample_seq.iter().map(|&c| c.to_ascii_uppercase()).collect())
                .unwrap();
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

    #[test]
    fn test_append_sample_with_paf_minus_strand_multi_segment() -> Result<()> {
        // ref = 5000 bp -> 2 segments (4096 + 904). Sample is RC(ref) with
        // SNPs in both segments plus an indel, exercised on '-' strand.
        let dir = tempfile::tempdir()?;
        let ref_seq = vec![b'A'; 5000];
        let ref_path = dir.path().join("ref.fa");
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", &ref_seq)]);

        // sample_fwd is almost identical to ref, with SNPs and one insertion.
        let mut sample_fwd = ref_seq.clone();
        sample_fwd[100] = b'C';
        sample_fwd[4100] = b'G';
        sample_fwd.insert(3000, b'T');

        let sample_seq: Vec<u8> = nt::rev_comp(&sample_fwd).collect();
        let sample_path = dir.path().join("sample.fa");
        write_fasta(sample_path.to_str().unwrap(), &[("chr1", &sample_seq)]);

        // PAF: '-' strand, CIGAR describes sample_fwd (== RC(sample)) vs ref.
        let paf_path = dir.path().join("sample.paf");
        write_paf(
            paf_path.to_str().unwrap(),
            &[paf_line(
                "chr1",
                5001,
                0,
                5001,
                "-",
                "chr1",
                5000,
                0,
                5000,
                "100=1X2899=1I1100=1X899=",
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
    fn test_encode_segment_lzdiff_empty_ref_groups_returns_error() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let ref_path = dir.path().join("ref.fa");
        write_fasta(ref_path.to_str().unwrap(), &[("chr1", b"ACGTACGTACGTACGT")]);
        let out_path = dir.path().join("out.pbit");

        let mut comp = Compressor::create(&out_path, ref_path.to_str().unwrap(), 4096, 15, 18)?;
        let seg = b"ACGT";
        let empty_ref_groups: &[u32] = &[];
        let result =
            comp.encode_segment_lzdiff("sample1", "chr1", 0, 0, seg, empty_ref_groups, false);

        assert!(result.is_err(), "expected error for empty ref_group_ids");
        let err = result.unwrap_err().to_string();
        assert!(
            err.contains("no reference groups"),
            "unexpected error: {}",
            err
        );
        Ok(())
    }
}
