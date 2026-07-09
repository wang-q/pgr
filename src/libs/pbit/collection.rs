//! Sample/contig/segment metadata (collection), serialized as fixed-size u32
//! LE fields + flate2 compression (no prefix coding).
//!
//! See `notes/design/pbit.md` §Sample Index for the binary layout.

use anyhow::{anyhow, Result};
use indexmap::IndexMap;

use super::format::{read_string, read_u32_le, write_string, write_u32_le};

/// Maximum sample count (1 million — far beyond any realistic archive).
const MAX_SAMPLE_COUNT: usize = 1_000_000;
/// Maximum contig count per sample (1 million).
const MAX_CONTIG_COUNT: usize = 1_000_000;
/// Maximum segment count per contig (100 million).
const MAX_SEGMENT_COUNT: usize = 100_000_000;

/// One segment of a sample's contig: a pointer into the reference group /
/// delta tables. `is_rev_comp` / `raw_length` / `encoding` live in
/// `DeltaEntry` (shared by all segments pointing to the same delta).
///
/// `ref_start` / `ref_end` are segment-relative offsets within the reference
/// 2bit record (used by CIGAR-encoded deltas; 0 for LZ-diff deltas).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SegmentDesc {
    pub ref_group_id: u32,
    pub delta_id: u32,
    pub ref_start: u32,
    pub ref_end: u32,
}

/// All segments of one contig within one sample.
#[derive(Debug, Clone)]
pub struct ContigSegs {
    pub contig_name: String,
    pub segments: Vec<SegmentDesc>,
}

/// Collection of samples and their contig segment mappings.
#[derive(Clone)]
pub struct Collection {
    /// sample name → contigs (each contig has a list of segments).
    pub samples: IndexMap<String, Vec<ContigSegs>>,
    pub cmd_line: String,
}

impl Collection {
    /// Create an empty collection.
    pub fn new() -> Self {
        Self {
            samples: IndexMap::new(),
            cmd_line: String::new(),
        }
    }

    /// Ensure `sample` is registered and has a contig entry for `contig`.
    /// Returns `false` if the sample already had this contig (duplicate).
    pub fn register_sample_contig(&mut self, sample: &str, contig: &str) -> bool {
        let contigs = self.samples.entry(sample.to_string()).or_default();
        if contigs.iter().any(|c| c.contig_name == contig) {
            return false;
        }
        contigs.push(ContigSegs {
            contig_name: contig.to_string(),
            segments: Vec::new(),
        });
        true
    }

    /// Ensure `sample` exists in the collection (with no contigs if new).
    pub fn ensure_sample(&mut self, sample: &str) {
        self.samples.entry(sample.to_string()).or_default();
    }

    /// Append a segment descriptor to `sample`'s `contig`.
    /// Panics-free: registers the sample/contig if missing.
    pub fn add_segment(
        &mut self,
        sample: &str,
        contig: &str,
        ref_group_id: u32,
        delta_id: u32,
        ref_start: u32,
        ref_end: u32,
    ) {
        self.register_sample_contig(sample, contig);
        let contigs = self.samples.get_mut(sample).expect("just registered");
        let entry = contigs
            .iter_mut()
            .find(|c| c.contig_name == contig)
            .expect("just registered");
        entry.segments.push(SegmentDesc {
            ref_group_id,
            delta_id,
            ref_start,
            ref_end,
        });
    }

    /// Return the segment list for `sample`'s `contig`, or `None`.
    pub fn get_contig_segments(&self, sample: &str, contig: &str) -> Option<&[SegmentDesc]> {
        self.samples
            .get(sample)?
            .iter()
            .find(|c| c.contig_name == contig)
            .map(|c| c.segments.as_slice())
    }

    /// List all sample names in insertion order.
    pub fn list_samples(&self) -> Vec<&str> {
        self.samples.keys().map(|s| s.as_str()).collect()
    }

    /// List all contig names for `sample` (insertion order), or empty if absent.
    pub fn list_contigs(&self, sample: &str) -> Vec<&str> {
        self.samples
            .get(sample)
            .map(|contigs| contigs.iter().map(|c| c.contig_name.as_str()).collect())
            .unwrap_or_default()
    }

    /// Number of registered samples.
    pub fn sample_count(&self) -> usize {
        self.samples.len()
    }

    /// Serialize to a flate2-compressed byte vector.
    ///
    /// Layout (before compression):
    /// ```text
    /// u32 sample_count
    /// for each sample:
    ///   u32 name_len + name_bytes
    ///   u32 contig_count
    ///   for each contig:
    ///     u32 contig_name_len + contig_name_bytes
    ///     u32 segment_count
    ///     for each segment:
    ///       u32 ref_group_id + u32 delta_id + u32 ref_start + u32 ref_end
    /// u32 cmd_line_len + cmd_line_bytes
    /// ```
    pub fn serialize(&self) -> Result<Vec<u8>> {
        let mut raw = Vec::new();
        write_u32_le(&mut raw, self.samples.len() as u32)?;
        for (sample, contigs) in &self.samples {
            write_string(&mut raw, sample)?;
            write_u32_le(&mut raw, contigs.len() as u32)?;
            for cs in contigs {
                write_string(&mut raw, &cs.contig_name)?;
                write_u32_le(&mut raw, cs.segments.len() as u32)?;
                for seg in &cs.segments {
                    write_u32_le(&mut raw, seg.ref_group_id)?;
                    write_u32_le(&mut raw, seg.delta_id)?;
                    write_u32_le(&mut raw, seg.ref_start)?;
                    write_u32_le(&mut raw, seg.ref_end)?;
                }
            }
        }
        write_string(&mut raw, &self.cmd_line)?;

        let mut encoder = flate2::write::GzEncoder::new(
            Vec::with_capacity(raw.len() / 2),
            flate2::Compression::default(),
        );
        use std::io::Write;
        encoder.write_all(&raw)?;
        Ok(encoder.finish()?)
    }

    /// Deserialize from a flate2-compressed byte vector (as produced by
    /// `serialize`).
    pub fn deserialize(data: &[u8]) -> Result<Self> {
        use std::io::Read;
        let mut decoder = flate2::read::GzDecoder::new(data);
        let mut raw = Vec::new();
        decoder.read_to_end(&mut raw)?;
        let mut cursor = std::io::Cursor::new(raw);

        let sample_count = read_u32_le(&mut cursor)? as usize;
        if sample_count > MAX_SAMPLE_COUNT {
            return Err(anyhow!(
                "sample_count {} exceeds maximum {}",
                sample_count,
                MAX_SAMPLE_COUNT
            ));
        }
        let mut samples: IndexMap<String, Vec<ContigSegs>> =
            IndexMap::with_capacity(sample_count.min(1024));
        for _ in 0..sample_count {
            let sample = read_string(&mut cursor)?;
            let contig_count = read_u32_le(&mut cursor)? as usize;
            if contig_count > MAX_CONTIG_COUNT {
                return Err(anyhow!(
                    "contig_count {} exceeds maximum {}",
                    contig_count,
                    MAX_CONTIG_COUNT
                ));
            }
            let mut contigs = Vec::with_capacity(contig_count.min(1024));
            for _ in 0..contig_count {
                let contig_name = read_string(&mut cursor)?;
                let segment_count = read_u32_le(&mut cursor)? as usize;
                if segment_count > MAX_SEGMENT_COUNT {
                    return Err(anyhow!(
                        "segment_count {} exceeds maximum {}",
                        segment_count,
                        MAX_SEGMENT_COUNT
                    ));
                }
                let mut segments = Vec::with_capacity(segment_count.min(1024));
                for _ in 0..segment_count {
                    let ref_group_id = read_u32_le(&mut cursor)?;
                    let delta_id = read_u32_le(&mut cursor)?;
                    let ref_start = read_u32_le(&mut cursor)?;
                    let ref_end = read_u32_le(&mut cursor)?;
                    segments.push(SegmentDesc {
                        ref_group_id,
                        delta_id,
                        ref_start,
                        ref_end,
                    });
                }
                contigs.push(ContigSegs {
                    contig_name,
                    segments,
                });
            }
            samples.insert(sample, contigs);
        }
        let cmd_line = read_string(&mut cursor)?;
        Ok(Self { samples, cmd_line })
    }
}

impl Default for Collection {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_collection_roundtrip() -> Result<()> {
        let col = Collection::new();
        let data = col.serialize()?;
        let back = Collection::deserialize(&data)?;
        assert_eq!(back.sample_count(), 0);
        assert!(back.cmd_line.is_empty());
        Ok(())
    }

    #[test]
    fn test_single_sample_single_contig() -> Result<()> {
        let mut col = Collection::new();
        col.add_segment("sample1", "chr1", 0, 0, 0, 0);
        col.add_segment("sample1", "chr1", 1, 2, 0, 0);

        let data = col.serialize()?;
        let back = Collection::deserialize(&data)?;
        assert_eq!(back.sample_count(), 1);
        assert_eq!(back.list_samples(), vec!["sample1"]);
        let segs = back.get_contig_segments("sample1", "chr1").unwrap();
        assert_eq!(segs.len(), 2);
        assert_eq!(
            segs[0],
            SegmentDesc {
                ref_group_id: 0,
                delta_id: 0,
                ref_start: 0,
                ref_end: 0
            }
        );
        assert_eq!(
            segs[1],
            SegmentDesc {
                ref_group_id: 1,
                delta_id: 2,
                ref_start: 0,
                ref_end: 0
            }
        );
        Ok(())
    }

    #[test]
    fn test_multiple_samples_contigs() -> Result<()> {
        let mut col = Collection::new();
        col.cmd_line = "pgr pbit create -r ref.fa".to_string();
        // sample1: chr1 (2 segments) + chr2 (1 segment)
        col.add_segment("sample1", "chr1", 0, 0, 0, 0);
        col.add_segment("sample1", "chr1", 1, 1, 0, 0);
        col.add_segment("sample1", "chr2", 2, 0, 0, 0);
        // sample2: chr1 (1 segment)
        col.add_segment("sample2", "chr1", 0, 3, 0, 0);

        let data = col.serialize()?;
        let back = Collection::deserialize(&data)?;
        assert_eq!(back.sample_count(), 2);
        assert_eq!(back.list_samples(), vec!["sample1", "sample2"]);
        assert_eq!(back.cmd_line, "pgr pbit create -r ref.fa");

        // sample1 contigs
        let s1_contigs = back.list_contigs("sample1");
        assert_eq!(s1_contigs, vec!["chr1", "chr2"]);
        let s1_chr1 = back.get_contig_segments("sample1", "chr1").unwrap();
        assert_eq!(s1_chr1.len(), 2);
        let s1_chr2 = back.get_contig_segments("sample1", "chr2").unwrap();
        assert_eq!(s1_chr2.len(), 1);
        assert_eq!(
            s1_chr2[0],
            SegmentDesc {
                ref_group_id: 2,
                delta_id: 0,
                ref_start: 0,
                ref_end: 0
            }
        );

        // sample2 contigs
        let s2_contigs = back.list_contigs("sample2");
        assert_eq!(s2_contigs, vec!["chr1"]);
        let s2_chr1 = back.get_contig_segments("sample2", "chr1").unwrap();
        assert_eq!(
            s2_chr1[0],
            SegmentDesc {
                ref_group_id: 0,
                delta_id: 3,
                ref_start: 0,
                ref_end: 0
            }
        );
        Ok(())
    }

    #[test]
    fn test_register_duplicate_contig() {
        let mut col = Collection::new();
        assert!(col.register_sample_contig("s1", "chr1"));
        assert!(!col.register_sample_contig("s1", "chr1")); // duplicate
        assert_eq!(col.list_contigs("s1").len(), 1);
    }

    #[test]
    fn test_get_missing_returns_none() {
        let col = Collection::new();
        assert!(col.get_contig_segments("nope", "chr1").is_none());
        let mut col2 = Collection::new();
        col2.add_segment("s1", "chr1", 0, 0, 0, 0);
        assert!(col2.get_contig_segments("s1", "chr2").is_none());
        assert!(col2.get_contig_segments("s2", "chr1").is_none());
    }

    #[test]
    fn test_unicode_sample_name() -> Result<()> {
        let mut col = Collection::new();
        col.add_segment("样本_1", "chr1", 0, 0, 0, 0);
        let data = col.serialize()?;
        let back = Collection::deserialize(&data)?;
        assert_eq!(back.list_samples(), vec!["样本_1"]);
        assert!(back.get_contig_segments("样本_1", "chr1").is_some());
        Ok(())
    }
}
