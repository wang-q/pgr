//! Compressor: builds a `.pbit` archive from reference + sample FASTA files.
//!
//! Holds a `W: Write + Seek` writer directly (no archive wrapper). The
//! reference layer is stored as standard 2bit records (reusing
//! `twobit::write_2bit_record`); sample segments are LZ-diff encoded against
//! the matching reference segment, flate2-compressed, and stored as delta
//! entries.

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::{Seek, SeekFrom, Write};
use std::path::Path;

use crate::libs::fmt::twobit::{read_2bit_record, write_2bit_record};
use crate::libs::nt;

use super::collection::Collection;
use super::decompressor::Decompressor;
use super::format::{
    read_u32_le, write_ref_index, write_u32_le, DeltaEntry, PbitFooter, PbitHeader, RefGroupEntry,
};
use super::segment::Segment;

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
    contig_ref_groups: HashMap<String, Vec<u32>>,
    /// Cached reference segment DNA (ASCII, forward) for rev-comp detection.
    /// Kept small: one segment_size-length string per ref_group.
    ref_seg_dna: Vec<Vec<u8>>,
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
            contig_ref_groups: HashMap::new(),
            ref_seg_dna: Vec::new(),
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
                let seg_str = std::str::from_utf8(seg).unwrap_or("");
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
                comp.ref_seg_dna.push(seg.to_vec());

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
        let mut ref_seg_dna: Vec<Vec<u8>> = Vec::with_capacity(ref_group_count);
        let mut contig_ref_groups: HashMap<String, Vec<u32>> = HashMap::new();
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
            ref_seg_dna.push(seq_bytes);
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
            ref_seg_dna,
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
            let ref_group_ids = match self.contig_ref_groups.get(contig_name) {
                Some(ids) => ids,
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
            let first_ref_dna = &self.ref_seg_dna[first_ref_group as usize];
            let is_rev_comp = detect_rev_comp(segs[0], first_ref_dna, self.kmer_len);

            for (seg_idx, seg) in segs.iter().enumerate() {
                // Match to reference segment by position (clamped to last).
                let ref_idx = seg_idx.min(ref_group_ids.len() - 1);
                let ref_group_id = ref_group_ids[ref_idx];

                // Apply rev-comp if detected.
                let encoded_seq: Vec<u8> = if is_rev_comp {
                    rev_comp_vec(seg)
                } else {
                    seg.to_vec()
                };

                // LZ-diff encode.
                let delta = self.segments[ref_group_id as usize].add(&encoded_seq);
                let raw_length = encoded_seq.len() as u32;

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
                        };
                        self.deltas[ref_group_id as usize].push(entry);
                        (self.deltas[ref_group_id as usize].len() - 1) as u32
                    }
                };

                self.collection
                    .add_segment(sample_name, contig_name, ref_group_id, delta_id);
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
fn flate2_compress(data: &[u8]) -> Result<Vec<u8>> {
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
}
