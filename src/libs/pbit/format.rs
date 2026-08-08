//! pbit file format constants, structures, and binary I/O.
//!
//! All integers use fixed-size little-endian encoding (no varint/prefix
//! coding). Strings use a u32 length prefix + UTF-8 bytes (no null
//! termination). See `notes/design/pbit.md` §文件格式规范 for the
//! complete specification.

use anyhow::{anyhow, Result};
use std::io::{Read, Seek, SeekFrom, Write};

/// pbit magic number: 'PBIT' in little-endian.
pub const PBIT_MAGIC: u32 = 0x54494250;
/// pbit format major version.
pub const PBIT_VERSION_MAJOR: u32 = 1;
/// pbit format minor version.
///
/// v5 (2026-08-05): sample `ContigSegs` gains `mask_blocks` (soft-mask
/// intervals, inherited from 2bit semantics), stored in the sample index.
///
/// v6 (2026-08-08): `DeltaEncoding::Raw` added; segments with no matching
/// reference content are stored verbatim so the archive is strictly lossless.
///
/// v7 (2026-08-08): `SegmentDesc` gains `q_start` and its `ref_start/ref_end`
/// become reference-file-global coordinates, so CIGAR deltas may reference any
/// reference interval (not just a single 4096 bp block). PAF-driven encoding
/// no longer requires segment-phase alignment.
///
/// v8 (2026-08-08): `DeltaEncoding::Raw` payload switches from
/// flate2(ASCII bases) to a standard 2bit record (same encoding as reference
/// segments; N via the 2bit N-block structure). 2bit is ~7 pp smaller than
/// flate2(ASCII) for DNA and needs no secondary compression.
///
/// v9 (2026-08-08): archives store per-sample PAF recovery data: verbatim
/// rows for records not CIGAR-encoded (small chains) plus the MAF score
/// (`ms:i`) of CIGAR-encoded records, so `pbit to-paf` can reproduce the
/// original PAF (big chains merged back to chain level, ms restored).
/// `SegmentDesc` gains `paf_id` (source PAF record id); footer gains
/// `paf_data_offset`.
///
/// v10 (2026-08-09): `DeltaEncoding::Identity` added; a segment exactly
/// identical to a reference interval (CIGAR ops all '=') is stored as a
/// zero-payload pointer to that interval instead of a packed CIGAR.
pub const PBIT_VERSION_MINOR: u32 = 10;
/// Current file version encoded as major*1000 + minor.
pub const PBIT_VERSION: u32 = PBIT_VERSION_MAJOR * 1000 + PBIT_VERSION_MINOR;

/// Delta encoding type stored in the on-disk delta header (10th byte).
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeltaEncoding {
    LzDiff = 0,
    Cigar = 1,
    /// Verbatim (flate2-compressed) sample sequence; used for segments with
    /// no matching reference content so the archive stays lossless.
    Raw = 2,
    /// Zero-payload pointer to a reference interval; used when the sample
    /// segment is exactly identical to it (CIGAR ops all '=').
    Identity = 3,
}

impl DeltaEncoding {
    pub fn from_u8(v: u8) -> Result<Self> {
        match v {
            0 => Ok(Self::LzDiff),
            1 => Ok(Self::Cigar),
            2 => Ok(Self::Raw),
            3 => Ok(Self::Identity),
            _ => Err(anyhow!("invalid DeltaEncoding: {}", v)),
        }
    }
}

/// PAF recovery data for one sample (v1009): (sample name, big-chain ms
/// table as (record_id, ms), verbatim small-chain PAF rows).
pub type PafRecovery = (String, Vec<(u32, i32)>, Vec<String>);

/// File header (fixed 36 bytes, at file start).
#[derive(Debug, Clone)]
pub struct PbitHeader {
    pub magic: u32,
    pub version: u32,
    pub segment_size: u32,
    pub kmer_len: u32,
    pub min_match_len: u32,
    pub ref_group_count: u32,
    pub sample_count: u32,
    pub ref_records_offset: u64,
}

impl PbitHeader {
    /// Create a new header with the given parameters and placeholder offsets.
    pub fn new(
        segment_size: u32,
        kmer_len: u32,
        min_match_len: u32,
        ref_group_count: u32,
        sample_count: u32,
    ) -> Self {
        Self {
            magic: PBIT_MAGIC,
            version: PBIT_VERSION,
            segment_size,
            kmer_len,
            min_match_len,
            ref_group_count,
            sample_count,
            ref_records_offset: 36,
        }
    }

    /// Read a 36-byte header from the current reader position.
    pub fn read_from<R: Read>(reader: &mut R) -> Result<Self> {
        let magic = read_u32_le(reader)?;
        if magic != PBIT_MAGIC {
            return Err(anyhow!(
                "Not a valid pbit file (magic: {:x}, expected {:x})",
                magic,
                PBIT_MAGIC
            ));
        }
        let version = read_u32_le(reader)?;
        if version != PBIT_VERSION {
            return Err(anyhow!(
                "unsupported pbit version {} (expected {})",
                version,
                PBIT_VERSION
            ));
        }
        let segment_size = read_u32_le(reader)?;
        let kmer_len = read_u32_le(reader)?;
        let min_match_len = read_u32_le(reader)?;
        let ref_group_count = read_u32_le(reader)?;
        let sample_count = read_u32_le(reader)?;
        let ref_records_offset = read_u64_le(reader)?;
        Ok(Self {
            magic,
            version,
            segment_size,
            kmer_len,
            min_match_len,
            ref_group_count,
            sample_count,
            ref_records_offset,
        })
    }

    /// Write the header (36 bytes) to the writer.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> Result<()> {
        writer.write_all(&self.magic.to_le_bytes())?;
        writer.write_all(&self.version.to_le_bytes())?;
        writer.write_all(&self.segment_size.to_le_bytes())?;
        writer.write_all(&self.kmer_len.to_le_bytes())?;
        writer.write_all(&self.min_match_len.to_le_bytes())?;
        writer.write_all(&self.ref_group_count.to_le_bytes())?;
        writer.write_all(&self.sample_count.to_le_bytes())?;
        writer.write_all(&self.ref_records_offset.to_le_bytes())?;
        Ok(())
    }
}

/// File footer (24 bytes; embedded reference index offsets live in the
/// Reference Index section's reference table).
///
/// Located by seeking to `file_size - 24`.
#[derive(Debug, Clone, Default)]
pub struct PbitFooter {
    pub ref_index_offset: u64,
    pub delta_data_offset: u64,
    pub sample_index_offset: u64,
    pub paf_data_offset: u64,
}

impl PbitFooter {
    /// Read a footer from the current reader position.
    pub fn read_from<R: Read>(reader: &mut R) -> Result<Self> {
        let ref_index_offset = read_u64_le(reader)?;
        let delta_data_offset = read_u64_le(reader)?;
        let sample_index_offset = read_u64_le(reader)?;
        let paf_data_offset = read_u64_le(reader)?;
        let footer = Self {
            ref_index_offset,
            delta_data_offset,
            sample_index_offset,
            paf_data_offset,
        };
        Ok(footer)
    }

    /// Write the footer (32 bytes) to the writer.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> Result<()> {
        writer.write_all(&self.ref_index_offset.to_le_bytes())?;
        writer.write_all(&self.delta_data_offset.to_le_bytes())?;
        writer.write_all(&self.sample_index_offset.to_le_bytes())?;
        writer.write_all(&self.paf_data_offset.to_le_bytes())?;
        Ok(())
    }

    /// Read the footer from the end of a seekable reader.
    pub fn read_at_end<R: Read + Seek>(reader: &mut R) -> Result<Self> {
        let file_size = reader.seek(SeekFrom::End(0))?;
        if file_size < 32 {
            return Err(anyhow!(
                "pbit file too small for footer: {} bytes",
                file_size
            ));
        }
        reader.seek(SeekFrom::Start(file_size - 32))?;
        Self::read_from(reader)
    }
}

/// Reference group index entry: one per reference segment (one 2bit record).
/// `ref_id` identifies the reference genome the segment belongs to.
#[derive(Debug, Clone)]
pub struct RefGroupEntry {
    pub contig_name: String,
    pub ref_id: u32,
    pub segment_offset: u64,
}

impl RefGroupEntry {
    /// Write a single reference group entry
    /// (u32 name_len + name + u32 ref_id + u64 offset).
    pub fn write_to<W: Write>(&self, writer: &mut W) -> Result<()> {
        write_string(writer, &self.contig_name)?;
        writer.write_all(&self.ref_id.to_le_bytes())?;
        writer.write_all(&self.segment_offset.to_le_bytes())?;
        Ok(())
    }

    /// Read a single reference group entry from the current reader position.
    pub fn read_from<R: Read>(reader: &mut R) -> Result<Self> {
        let contig_name = read_string(reader)?;
        let ref_id = read_u32_le(reader)?;
        let segment_offset = read_u64_le(reader)?;
        Ok(Self {
            contig_name,
            ref_id,
            segment_offset,
        })
    }
}

/// Write a reference index section: u32 ref_group_count + entries.
pub fn write_ref_index<W: Write>(writer: &mut W, entries: &[RefGroupEntry]) -> Result<()> {
    writer.write_all(&(entries.len() as u32).to_le_bytes())?;
    for entry in entries {
        entry.write_to(writer)?;
    }
    Ok(())
}

/// Read a reference index section and return the entries.
pub fn read_ref_index<R: Read>(reader: &mut R) -> Result<Vec<RefGroupEntry>> {
    let count = read_u32_le(reader)? as usize;
    if count > MAX_REF_GROUPS {
        return Err(anyhow!(
            "reference index count {} exceeds maximum {}; corrupt or malicious",
            count,
            MAX_REF_GROUPS
        ));
    }
    let mut entries = Vec::with_capacity(count);
    for _ in 0..count {
        entries.push(RefGroupEntry::read_from(reader)?);
    }
    Ok(entries)
}

/// One reference genome in the archive: its name and segment group range.
#[derive(Debug, Clone, Default)]
pub struct RefTableEntry {
    pub ref_name: String,
    pub group_start: u32,
    pub group_count: u32,
}

/// Write the reference table: u32 ref_count + entries (appended after the
/// per-group entries in the Reference Index section).
pub fn write_ref_table<W: Write>(writer: &mut W, refs: &[RefTableEntry]) -> Result<()> {
    writer.write_all(&(refs.len() as u32).to_le_bytes())?;
    for r in refs {
        write_string(writer, &r.ref_name)?;
        writer.write_all(&r.group_start.to_le_bytes())?;
        writer.write_all(&r.group_count.to_le_bytes())?;
    }
    Ok(())
}

/// Read a reference table (after the per-group entries).
pub fn read_ref_table<R: Read>(reader: &mut R) -> Result<Vec<RefTableEntry>> {
    let count = read_u32_le(reader)? as usize;
    if count > MAX_REF_GROUPS {
        return Err(anyhow!(
            "reference table count {} exceeds maximum {}; corrupt or malicious",
            count,
            MAX_REF_GROUPS
        ));
    }
    let mut refs = Vec::with_capacity(count);
    for _ in 0..count {
        let ref_name = read_string(reader)?;
        let group_start = read_u32_le(reader)?;
        let group_count = read_u32_le(reader)?;
        refs.push(RefTableEntry {
            ref_name,
            group_start,
            group_count,
        });
    }
    Ok(refs)
}

/// In-memory delta header (10 bytes): properties of a delta encoding, loaded
/// during `Decompressor::new` without the packed data.
#[derive(Debug, Clone, Copy)]
pub struct DeltaMeta {
    pub is_rev_comp: bool,
    pub raw_length: u32,
    pub packed_size: u32,
    pub encoding: DeltaEncoding,
}

impl DeltaMeta {
    /// Write the 10-byte delta header (no packed data).
    pub fn write_header<W: Write>(&self, writer: &mut W) -> Result<()> {
        writer.write_all(&[self.is_rev_comp as u8])?;
        writer.write_all(&self.raw_length.to_le_bytes())?;
        writer.write_all(&self.packed_size.to_le_bytes())?;
        writer.write_all(&[self.encoding as u8])?;
        Ok(())
    }

    /// Read the 10-byte delta header (no packed data).
    pub fn read_header<R: Read>(reader: &mut R) -> Result<Self> {
        let mut byte = [0u8; 1];
        reader.read_exact(&mut byte)?;
        let is_rev_comp = byte[0] != 0;
        let raw_length = read_u32_le(reader)?;
        let packed_size = read_u32_le(reader)?;
        reader.read_exact(&mut byte)?;
        let encoding = DeltaEncoding::from_u8(byte[0])?;
        Ok(Self {
            is_rev_comp,
            raw_length,
            packed_size,
            encoding,
        })
    }
}

/// Full delta entry: header + packed (flate2-compressed) data.
#[derive(Debug, Clone)]
pub struct DeltaEntry {
    pub is_rev_comp: bool,
    pub raw_length: u32,
    pub packed_data: Vec<u8>,
    pub encoding: DeltaEncoding,
}

impl DeltaEntry {
    /// Return the in-memory metadata header for this delta entry.
    pub fn meta(&self) -> DeltaMeta {
        DeltaMeta {
            is_rev_comp: self.is_rev_comp,
            raw_length: self.raw_length,
            packed_size: self.packed_data.len() as u32,
            encoding: self.encoding,
        }
    }

    /// Write the full delta entry (10-byte header + packed data).
    pub fn write_to<W: Write>(&self, writer: &mut W) -> Result<()> {
        self.meta().write_header(writer)?;
        writer.write_all(&self.packed_data)?;
        Ok(())
    }

    /// Read the full delta entry (10-byte header + packed data) from the
    /// current reader position.
    pub fn read_from<R: Read>(reader: &mut R) -> Result<Self> {
        let meta = DeltaMeta::read_header(reader)?;
        if meta.packed_size as usize > MAX_PACKED_SIZE {
            return Err(anyhow!(
                "delta packed_size {} exceeds maximum {}; corrupt or malicious",
                meta.packed_size,
                MAX_PACKED_SIZE
            ));
        }
        let mut packed_data = vec![0u8; meta.packed_size as usize];
        reader.read_exact(&mut packed_data)?;
        Ok(Self {
            is_rev_comp: meta.is_rev_comp,
            raw_length: meta.raw_length,
            packed_data,
            encoding: meta.encoding,
        })
    }
}

/// Write a u32 little-endian.
pub fn write_u32_le<W: Write>(writer: &mut W, val: u32) -> Result<()> {
    writer.write_all(&val.to_le_bytes())?;
    Ok(())
}

/// Write a u64 little-endian.
pub fn write_u64_le<W: Write>(writer: &mut W, val: u64) -> Result<()> {
    writer.write_all(&val.to_le_bytes())?;
    Ok(())
}

/// Read a u32 little-endian.
pub fn read_u32_le<R: Read>(reader: &mut R) -> Result<u32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

/// Read a u64 little-endian.
pub fn read_u64_le<R: Read>(reader: &mut R) -> Result<u64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

/// Write a length-prefixed string (u32 len + UTF-8 bytes).
pub fn write_string<W: Write>(writer: &mut W, s: &str) -> Result<()> {
    let bytes = s.as_bytes();
    writer.write_all(&(bytes.len() as u32).to_le_bytes())?;
    writer.write_all(bytes)?;
    Ok(())
}

/// Maximum string length (16 MB) — guards against malicious length prefixes.
const MAX_STRING_LEN: usize = 16 * 1024 * 1024;

/// Maximum reference groups (segments) read from an index or header. A
/// malicious count drives `Vec::with_capacity` before validation, so it must
/// be bounded to avoid a multi-GB upfront allocation. 1M segments covers
/// any realistic reference (whole human genome ≈ 750k segments at 4096 bp).
pub(crate) const MAX_REF_GROUPS: usize = 1_000_000;

/// Maximum per-group delta count. Bounds the per-group `Vec::with_capacity`
/// in the decompressor's delta scan.
pub(crate) const MAX_DELTAS_PER_GROUP: usize = 1_000_000;

/// Maximum packed (compressed) size of a single delta. A delta encodes at
/// most one reference segment, whose compressed form is far below this bound
/// for any realistic segment size (default 4096 bp). An inflated value read
/// from a malicious archive would otherwise trigger a multi-GB allocation in
/// `Decompressor::decode_delta` / `DeltaEntry::read_from`.
pub const MAX_PACKED_SIZE: usize = 256 * 1024 * 1024;

/// Maximum uncompressed size of a single delta's payload (LZ-diff raw delta
/// or packed CIGAR). A delta decodes at most one segment, so its uncompressed
/// form is O(segment_size) and far below this bound for any realistic archive.
/// An attacker can embed a gzip bomb (small compressed, huge expansion) inside
/// an otherwise small delta; bounding the decompressed size prevents a
/// multi-GB allocation in `Decompressor::decode_delta` / `unpack_cigar`.
pub(crate) const MAX_DELTA_UNCOMPRESSED: usize = 256 * 1024 * 1024;

/// Read a length-prefixed string (u32 len + UTF-8 bytes).
pub fn read_string<R: Read>(reader: &mut R) -> Result<String> {
    let len = read_u32_le(reader)? as usize;
    if len > MAX_STRING_LEN {
        return Err(anyhow!(
            "string length {} exceeds maximum {}",
            len,
            MAX_STRING_LEN
        ));
    }
    let mut buf = vec![0u8; len];
    reader.read_exact(&mut buf)?;
    Ok(String::from_utf8(buf)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_header_roundtrip() -> Result<()> {
        let header = PbitHeader::new(4096, 15, 18, 100, 10);
        let mut buf = Vec::new();
        header.write_to(&mut buf)?;
        assert_eq!(buf.len(), 36);

        let mut cursor = Cursor::new(buf);
        let read = PbitHeader::read_from(&mut cursor)?;
        assert_eq!(read.magic, PBIT_MAGIC);
        assert_eq!(read.version, PBIT_VERSION);
        assert_eq!(read.segment_size, 4096);
        assert_eq!(read.kmer_len, 15);
        assert_eq!(read.min_match_len, 18);
        assert_eq!(read.ref_group_count, 100);
        assert_eq!(read.sample_count, 10);
        assert_eq!(read.ref_records_offset, 36);
        Ok(())
    }

    #[test]
    fn test_header_bad_magic() {
        let mut buf = Vec::new();
        buf.extend_from_slice(&0xDEADBEEFu32.to_le_bytes());
        buf.extend_from_slice(&[0u8; 32]);

        let mut cursor = Cursor::new(buf);
        let res = PbitHeader::read_from(&mut cursor);
        assert!(res.is_err());
        assert!(res
            .unwrap_err()
            .to_string()
            .contains("Not a valid pbit file"));
    }

    #[test]
    fn test_header_bad_version() {
        let mut buf = Vec::new();
        buf.extend_from_slice(&PBIT_MAGIC.to_le_bytes());
        buf.extend_from_slice(&(PBIT_VERSION + 1).to_le_bytes());
        buf.extend_from_slice(&[0u8; 28]);

        let mut cursor = Cursor::new(buf);
        let res = PbitHeader::read_from(&mut cursor);
        assert!(res.is_err());
        assert!(res
            .unwrap_err()
            .to_string()
            .contains("unsupported pbit version"));
    }

    #[test]
    fn test_footer_roundtrip() -> Result<()> {
        let footer = PbitFooter {
            ref_index_offset: 1024,
            delta_data_offset: 2048,
            sample_index_offset: 4096,
            paf_data_offset: 0,
        };
        let mut buf = Vec::new();
        footer.write_to(&mut buf)?;
        assert_eq!(buf.len(), 32);

        let mut cursor = Cursor::new(buf);
        let read = PbitFooter::read_from(&mut cursor)?;
        assert_eq!(read.ref_index_offset, 1024);
        assert_eq!(read.delta_data_offset, 2048);
        assert_eq!(read.sample_index_offset, 4096);
        assert_eq!(read.paf_data_offset, 0);
        Ok(())
    }

    #[test]
    fn test_footer_read_at_end() -> Result<()> {
        let footer = PbitFooter {
            ref_index_offset: 100,
            delta_data_offset: 200,
            sample_index_offset: 300,
            paf_data_offset: 0,
        };
        // Prepend some padding before the footer
        let mut buf = vec![0xABu8; 50];
        let mut footer_buf = Vec::new();
        footer.write_to(&mut footer_buf)?;
        buf.extend_from_slice(&footer_buf);

        let mut cursor = Cursor::new(buf);
        let read = PbitFooter::read_at_end(&mut cursor)?;
        assert_eq!(read.ref_index_offset, 100);
        assert_eq!(read.delta_data_offset, 200);
        assert_eq!(read.sample_index_offset, 300);
        Ok(())
    }

    #[test]
    fn test_ref_group_entry_roundtrip() -> Result<()> {
        let entry = RefGroupEntry {
            contig_name: "chr1".to_string(),
            ref_id: 0,
            segment_offset: 12345,
        };
        let mut buf = Vec::new();
        entry.write_to(&mut buf)?;

        let mut cursor = Cursor::new(buf);
        let read = RefGroupEntry::read_from(&mut cursor)?;
        assert_eq!(read.contig_name, "chr1");
        assert_eq!(read.segment_offset, 12345);
        Ok(())
    }

    #[test]
    fn test_ref_index_roundtrip() -> Result<()> {
        let entries = vec![
            RefGroupEntry {
                contig_name: "chr1".to_string(),
                ref_id: 0,
                segment_offset: 100,
            },
            RefGroupEntry {
                contig_name: "chr2".to_string(),
                ref_id: 0,
                segment_offset: 200,
            },
        ];
        let mut buf = Vec::new();
        write_ref_index(&mut buf, &entries)?;

        let mut cursor = Cursor::new(buf);
        let read = read_ref_index(&mut cursor)?;
        assert_eq!(read.len(), 2);
        assert_eq!(read[0].contig_name, "chr1");
        assert_eq!(read[0].segment_offset, 100);
        assert_eq!(read[1].contig_name, "chr2");
        assert_eq!(read[1].segment_offset, 200);
        Ok(())
    }

    #[test]
    fn test_delta_meta_roundtrip() -> Result<()> {
        let meta = DeltaMeta {
            is_rev_comp: true,
            raw_length: 4096,
            packed_size: 512,
            encoding: DeltaEncoding::LzDiff,
        };
        let mut buf = Vec::new();
        meta.write_header(&mut buf)?;
        assert_eq!(buf.len(), 10);

        let mut cursor = Cursor::new(buf);
        let read = DeltaMeta::read_header(&mut cursor)?;
        assert!(read.is_rev_comp);
        assert_eq!(read.raw_length, 4096);
        assert_eq!(read.packed_size, 512);
        assert_eq!(read.encoding, DeltaEncoding::LzDiff);
        Ok(())
    }

    #[test]
    fn test_delta_entry_roundtrip() -> Result<()> {
        let entry = DeltaEntry {
            is_rev_comp: false,
            raw_length: 100,
            packed_data: vec![1, 2, 3, 4, 5],
            encoding: DeltaEncoding::Cigar,
        };
        let mut buf = Vec::new();
        entry.write_to(&mut buf)?;

        let mut cursor = Cursor::new(buf);
        let read = DeltaEntry::read_from(&mut cursor)?;
        assert!(!read.is_rev_comp);
        assert_eq!(read.raw_length, 100);
        assert_eq!(read.packed_data, vec![1, 2, 3, 4, 5]);
        assert_eq!(read.encoding, DeltaEncoding::Cigar);
        Ok(())
    }

    #[test]
    fn test_delta_entry_identity_roundtrip() -> Result<()> {
        // Identity entries carry no payload; the 10-byte header alone must
        // roundtrip with encoding preserved.
        let entry = DeltaEntry {
            is_rev_comp: true,
            raw_length: 4096,
            packed_data: Vec::new(),
            encoding: DeltaEncoding::Identity,
        };
        let mut buf = Vec::new();
        entry.write_to(&mut buf)?;
        assert_eq!(buf.len(), 10);

        let mut cursor = Cursor::new(buf);
        let read = DeltaEntry::read_from(&mut cursor)?;
        assert!(read.is_rev_comp);
        assert_eq!(read.raw_length, 4096);
        assert!(read.packed_data.is_empty());
        assert_eq!(read.encoding, DeltaEncoding::Identity);
        Ok(())
    }

    #[test]
    fn test_string_roundtrip() -> Result<()> {
        let s = "hello pbit 世界";
        let mut buf = Vec::new();
        write_string(&mut buf, s)?;

        let mut cursor = Cursor::new(buf);
        let read = read_string(&mut cursor)?;
        assert_eq!(read, s);
        Ok(())
    }

    #[test]
    fn test_empty_pbit_roundtrip() -> Result<()> {
        // An empty .pbit file: Header (placeholder offsets) + Footer (zero offsets).
        // No reference records, no index, no delta data, no sample index.
        let header = PbitHeader::new(4096, 15, 18, 0, 0);
        let footer = PbitFooter {
            ref_index_offset: 0,
            delta_data_offset: 0,
            sample_index_offset: 0,
            paf_data_offset: 0,
        };

        let mut buf = Vec::new();
        header.write_to(&mut buf)?;
        // No reference records (ref_group_count = 0)
        // No ref index / delta data / sample index sections
        footer.write_to(&mut buf)?;

        assert_eq!(buf.len(), 68); // 36 + 32

        let mut cursor = Cursor::new(buf);
        let read_header = PbitHeader::read_from(&mut cursor)?;
        assert_eq!(read_header.ref_group_count, 0);
        assert_eq!(read_header.sample_count, 0);

        let read_footer = PbitFooter::read_at_end(&mut cursor)?;
        assert_eq!(read_footer.ref_index_offset, 0);
        assert_eq!(read_footer.delta_data_offset, 0);
        assert_eq!(read_footer.sample_index_offset, 0);
        Ok(())
    }

    #[test]
    fn test_minimal_pbit_roundtrip() -> Result<()> {
        // A minimal .pbit: Header + empty ref index + empty delta data + empty
        // sample index + Footer.
        use crate::libs::fmt::twobit::write_2bit_record;

        let header = PbitHeader::new(4096, 15, 18, 1, 0);
        let mut buf = Vec::new();

        // Header (36 bytes)
        header.write_to(&mut buf)?;

        // Reference Records: one 2bit record (segment of "chr1")
        let ref_offset = buf.len() as u64;
        write_2bit_record(&mut buf, "ACGTACGT", false)?;

        // Reference Index
        let ref_index_offset = buf.len() as u64;
        let entries = vec![RefGroupEntry {
            contig_name: "chr1".to_string(),
            ref_id: 0,
            segment_offset: ref_offset,
        }];
        write_ref_index(&mut buf, &entries)?;

        // Delta Data: ref_group_count=1, delta_count=0
        let delta_data_offset = buf.len() as u64;
        buf.extend_from_slice(&1u32.to_le_bytes()); // ref_group_count
        buf.extend_from_slice(&0u32.to_le_bytes()); // delta_count

        // Sample Index: sample_count=0, cmd_line_len=0
        let sample_index_offset = buf.len() as u64;
        buf.extend_from_slice(&0u32.to_le_bytes()); // sample_count
        buf.extend_from_slice(&0u32.to_le_bytes()); // cmd_line_len

        // Footer
        let footer = PbitFooter {
            ref_index_offset,
            delta_data_offset,
            sample_index_offset,
            paf_data_offset: 0,
        };
        footer.write_to(&mut buf)?;

        // Read back
        let mut cursor = Cursor::new(buf);
        let read_header = PbitHeader::read_from(&mut cursor)?;
        assert_eq!(read_header.ref_group_count, 1);
        assert_eq!(read_header.sample_count, 0);

        let read_footer = PbitFooter::read_at_end(&mut cursor)?;

        // Read reference index
        cursor.seek(SeekFrom::Start(read_footer.ref_index_offset))?;
        let ref_entries = read_ref_index(&mut cursor)?;
        assert_eq!(ref_entries.len(), 1);
        assert_eq!(ref_entries[0].contig_name, "chr1");

        // Read delta data header
        cursor.seek(SeekFrom::Start(read_footer.delta_data_offset))?;
        let ref_group_count = read_u32_le(&mut cursor)?;
        assert_eq!(ref_group_count, 1);
        let delta_count = read_u32_le(&mut cursor)?;
        assert_eq!(delta_count, 0);

        // Read sample index
        cursor.seek(SeekFrom::Start(read_footer.sample_index_offset))?;
        let sample_count = read_u32_le(&mut cursor)?;
        assert_eq!(sample_count, 0);

        Ok(())
    }
}
