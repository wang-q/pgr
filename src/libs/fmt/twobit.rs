use anyhow::{anyhow, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom, Write};
use std::ops::{Bound, Deref, Range, RangeBounds};
use std::path::Path;

const TWOBIT_MAGIC: u32 = 0x1A412743;
const TWOBIT_MAGIC_SWAPPED: u32 = 0x4327411A;

/// A half-open interval representing a masked region within a sequence.
pub type Block = Range<usize>;

/// A collection of masked regions (N-blocks or soft-mask blocks) for a sequence.
#[derive(Debug, Clone)]
pub struct Blocks(pub Vec<Block>);

impl Deref for Blocks {
    type Target = [Block];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl Blocks {
    /// Return an iterator over blocks overlapping the given half-open range.
    pub fn overlaps(&self, range: impl RangeBounds<usize>) -> impl Iterator<Item = &Block> {
        let start = match range.start_bound() {
            Bound::Included(&x) => x,
            Bound::Excluded(&x) => x + 1,
            Bound::Unbounded => 0,
        };

        let end = match range.end_bound() {
            Bound::Included(&x) => Bound::Included(x),
            Bound::Excluded(&x) => Bound::Excluded(x),
            Bound::Unbounded => Bound::Unbounded,
        };

        self.0
            .iter()
            // Assume blocks are sorted
            .skip_while(move |block| block.end <= start)
            .take_while(move |block| match end {
                Bound::Included(e) => block.start <= e,
                Bound::Excluded(e) => block.start < e,
                Bound::Unbounded => true,
            })
    }

    /// Replace overlapping positions in `seq` with 'N', using 0-based global coordinates.
    pub fn apply_hard_mask(&self, seq: &mut [u8], offset: usize) {
        let (start, end) = (offset, offset + seq.len());
        for block in self.overlaps(start..end) {
            let seq_start = start.max(block.start) - start;
            let seq_end = end.min(block.end) - start;
            seq[seq_start..seq_end].fill(b'N');
        }
    }

    /// Convert overlapping positions in `seq` to lowercase, using 0-based global coordinates.
    pub fn apply_soft_mask(&self, seq: &mut [u8], offset: usize) {
        let (start, end) = (offset, offset + seq.len());
        for block in self.overlaps(start..end) {
            let seq_start = start.max(block.start) - start;
            let seq_end = end.min(block.end) - start;
            for byte in &mut seq[seq_start..seq_end] {
                *byte = byte.to_ascii_lowercase();
            }
        }
    }

    /// Pack a DNA string into 2-bit data and extract N-blocks and soft-mask blocks.
    ///
    /// Returns `(packed_dna, n_blocks, mask_blocks, dna_size)`.
    pub fn from_dna(dna: &str, do_mask: bool) -> (Vec<u8>, Blocks, Blocks, u32) {
        let len = dna.len();
        let mut n_blocks_vec = Vec::new();
        let mut mask_blocks_vec = Vec::new();
        let mut packed_dna = Vec::with_capacity(len.div_ceil(4));

        // Temporary buffer for current byte construction
        let mut current_byte = 0u8;
        let mut bit_offset = 6; // Starts at 6 (highest bits), then 4, 2, 0

        // State tracking for blocks
        let mut in_n = false;
        let mut n_start = 0;
        let mut in_mask = false;
        let mut mask_start = 0;

        for (i, c) in dna.chars().enumerate() {
            // Handle N-blocks (Hard mask)
            // Any character that is not a valid DNA base (A/C/G/T) is treated as N,
            // matching UCSC faToTwoBit semantics.
            let is_n = !matches!(c.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T');
            if is_n {
                if !in_n {
                    in_n = true;
                    n_start = i;
                }
            } else if in_n {
                in_n = false;
                n_blocks_vec.push(n_start..i);
            }

            // Handle Mask-blocks (Soft mask - lowercase)
            // Note: Ns are usually not counted as soft-mask in UCSC,
            // but if they are lowercase 'n', they might be?
            // twoBit.c: "lower-case characters are masked".
            // faToTwoBit.c: unknownToN converts to 'N' or 'n' based on case.
            // But usually N-mask takes precedence or is separate.
            // Let's stick to simple logic: if it's lowercase, it's a mask block.
            // However, typical usage is: valid bases in lowercase -> mask.
            let is_lower = c.is_ascii_lowercase();
            if do_mask && is_lower {
                if !in_mask {
                    in_mask = true;
                    mask_start = i;
                }
            } else if in_mask {
                in_mask = false;
                mask_blocks_vec.push(mask_start..i);
            }

            // Pack DNA
            // T=00, C=01, A=10, G=11
            // N is treated as T (00) usually, or C?
            // UCSC twoBitFromDnaSeq: val = ntVal[c]
            // We need a map. T/t/N/n -> ?
            // Usually arbitrary for N since it's masked. Let's use T (00).
            let val = match c {
                'T' | 't' => 0,
                'C' | 'c' => 1,
                'A' | 'a' => 2,
                'G' | 'g' => 3,
                _ => 0, // Treat N and others as T
            };

            current_byte |= val << bit_offset;

            if bit_offset == 0 {
                packed_dna.push(current_byte);
                current_byte = 0;
                bit_offset = 6;
            } else {
                bit_offset -= 2;
            }
        }

        // Close open blocks
        if in_n {
            n_blocks_vec.push(n_start..len);
        }
        if in_mask {
            mask_blocks_vec.push(mask_start..len);
        }

        // Push last partial byte if exists
        if bit_offset != 6 {
            packed_dna.push(current_byte);
        }

        (
            packed_dna,
            Blocks(n_blocks_vec),
            Blocks(mask_blocks_vec),
            len as u32,
        )
    }
}

/// Read a single 2bit record from the current reader position and return the
/// decoded DNA string with masks applied. Reused by `TwoBitFile` and `pbit::Decompressor`.
/// The reader must be positioned at the start of the record (dna_size field).
pub fn read_2bit_record<R: Read + Seek>(
    reader: &mut R,
    is_swapped: bool,
    start: Option<usize>,
    end: Option<usize>,
    no_mask: bool,
) -> Result<String> {
    let dna_size = read_u32(reader, is_swapped)? as usize;

    let n_blocks = read_blocks(reader, is_swapped)?;
    let mask_blocks = read_blocks(reader, is_swapped)?;

    let _reserved = read_u32(reader, is_swapped)?;

    let start_pos = start.unwrap_or(0);
    let end_pos = end.unwrap_or(dna_size).min(dna_size);

    if start_pos >= end_pos {
        return Ok(String::new());
    }

    // Currently at the start of packed DNA.
    let packed_dna_start = reader.stream_position()?;

    let first_byte_idx = start_pos / 4;
    let last_byte_idx = (end_pos - 1) / 4;

    reader.seek(SeekFrom::Start(packed_dna_start + first_byte_idx as u64))?;

    let mut packed_buf = vec![0u8; last_byte_idx - first_byte_idx + 1];
    reader.read_exact(&mut packed_buf)?;

    let table = [b'T', b'C', b'A', b'G'];

    let mut seq_vec = Vec::with_capacity(end_pos - start_pos);
    for i in start_pos..end_pos {
        let global_byte_idx = i / 4;
        let local_byte_idx = global_byte_idx - first_byte_idx;
        let bit_offset = 6 - 2 * (i % 4); // 0->6, 1->4, 2->2, 3->0

        let byte = packed_buf[local_byte_idx];
        let val = (byte >> bit_offset) & 3;
        seq_vec.push(table[val as usize]);
    }

    n_blocks.apply_hard_mask(&mut seq_vec, start_pos);

    if !no_mask {
        mask_blocks.apply_soft_mask(&mut seq_vec, start_pos);
    }

    String::from_utf8(seq_vec).map_err(|e| anyhow!("invalid UTF-8 in 2bit sequence: {}", e))
}

/// Write a single 2bit record (dna_size + n_blocks + mask_blocks + reserved +
/// packed_dna) to the writer. Reused by `TwoBitWriter` and `pbit::Compressor`.
pub fn write_2bit_record<W: Write>(writer: &mut W, dna: &str, do_mask: bool) -> Result<()> {
    let (packed, n_blocks, mask_blocks, dna_size) = Blocks::from_dna(dna, do_mask);

    writer.write_all(&dna_size.to_le_bytes())?;

    // Write N Blocks
    writer.write_all(&(n_blocks.0.len() as u32).to_le_bytes())?;
    for block in &n_blocks.0 {
        writer.write_all(&(block.start as u32).to_le_bytes())?;
    }
    for block in &n_blocks.0 {
        writer.write_all(&((block.end - block.start) as u32).to_le_bytes())?;
    }

    // Write Mask Blocks
    writer.write_all(&(mask_blocks.0.len() as u32).to_le_bytes())?;
    for block in &mask_blocks.0 {
        writer.write_all(&(block.start as u32).to_le_bytes())?;
    }
    for block in &mask_blocks.0 {
        writer.write_all(&((block.end - block.start) as u32).to_le_bytes())?;
    }

    // Reserved
    writer.write_all(&0u32.to_le_bytes())?;

    // Packed DNA
    writer.write_all(&packed)?;

    Ok(())
}

/// Read a block list (count + starts + sizes) from a 2bit record.
pub fn read_blocks<R: Read>(reader: &mut R, is_swapped: bool) -> Result<Blocks> {
    let count = read_u32(reader, is_swapped)? as usize;
    let starts = read_u32_vec(reader, count, is_swapped)?;
    let sizes = read_u32_vec(reader, count, is_swapped)?;

    let blocks: Vec<Block> = starts
        .into_iter()
        .zip(sizes)
        .map(|(start, size)| {
            let s = start as usize;
            let e = s + size as usize;
            s..e
        })
        .collect();

    Ok(Blocks(blocks))
}

/// A writer for the 2bit binary sequence format.
pub struct TwoBitWriter<W> {
    writer: W,
}

impl<W: std::io::Write> TwoBitWriter<W> {
    /// Create a new 2bit writer wrapping the given output.
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    /// Write a complete 2bit file containing the given sequences.
    ///
    /// `do_mask` controls whether lowercase bases are stored as soft-mask blocks.
    pub fn write(&mut self, sequences: &[(&str, &str)], do_mask: bool) -> Result<()> {
        // 1. Write Header
        self.writer.write_all(&TWOBIT_MAGIC.to_le_bytes())?;
        self.writer.write_all(&1u32.to_le_bytes())?; // Version 1
        self.writer
            .write_all(&(sequences.len() as u32).to_le_bytes())?; // SeqCount
        self.writer.write_all(&0u32.to_le_bytes())?; // Reserved

        // 2. Calculate offsets and Write Index
        // Header is 16 bytes; index starts at byte 16.
        // Index entry: NameLen(1) + Name(N) + Offset(8)

        let mut current_offset = 16u64;
        for (name, _) in sequences {
            current_offset += 1 + name.len() as u64 + 8;
        }
        // current_offset is now the offset of the first sequence record.
        let mut record_offsets = Vec::new();
        let mut running_offset = current_offset;

        for (name, dna) in sequences {
            // Write Index Entry
            let name_bytes = name.as_bytes();
            if name_bytes.len() > 255 {
                return Err(anyhow!("Sequence name too long: {}", name));
            }
            self.writer.write_all(&[name_bytes.len() as u8])?;
            self.writer.write_all(name_bytes)?;
            self.writer.write_all(&running_offset.to_le_bytes())?;

            record_offsets.push(running_offset);

            // Calculate next offset
            // Record overhead:
            // size(4) + nBlockCount(4) + nStarts(...) + nSizes(...) +
            // maskBlockCount(4) + maskStarts(...) + maskSizes(...) +
            // reserved(4) + packedDna(...)

            let (_, n_blocks, mask_blocks, _) = Blocks::from_dna(dna, do_mask);

            let n_count = n_blocks.0.len() as u64;
            let mask_count = mask_blocks.0.len() as u64;
            let packed_len = dna.len().div_ceil(4) as u64;

            let record_size = 4 + // dnaSize
                4 + (n_count * 4) + (n_count * 4) + // N blocks
                4 + (mask_count * 4) + (mask_count * 4) + // Mask blocks
                4 + // reserved
                packed_len;

            running_offset += record_size;
        }

        // 3. Write Records
        for (_, dna) in sequences.iter() {
            write_2bit_record(&mut self.writer, dna, do_mask)?;
        }

        Ok(())
    }
}

/// A random-access reader for the 2bit binary sequence format.
#[derive(Debug)]
pub struct TwoBitFile<R> {
    reader: R,
    pub sequence_offsets: HashMap<String, u64>,
    pub sequence_order: Vec<String>,
    is_swapped: bool,
    pub version: u32,
}

impl TwoBitFile<BufReader<File>> {
    /// Open a 2bit file from a path using a buffered file reader.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        Self::new(reader)
    }
}

impl TwoBitFile<std::io::Cursor<Vec<u8>>> {
    /// Read an entire 2bit file into memory and open it from a cursor.
    pub fn open_and_read<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mut buf = vec![];
        File::open(path)?.read_to_end(&mut buf)?;
        Self::new(std::io::Cursor::new(buf))
    }
}

impl<R: Read + Seek> TwoBitFile<R> {
    /// Parse a 2bit header and index from the given reader.
    pub fn new(mut reader: R) -> Result<Self> {
        // Read magic
        let mut buf = [0u8; 4];
        reader.read_exact(&mut buf)?;
        let magic = u32::from_le_bytes(buf);

        let is_swapped = if magic == TWOBIT_MAGIC {
            false
        } else if magic == TWOBIT_MAGIC_SWAPPED {
            true
        } else {
            return Err(anyhow!("Not a valid 2bit file (magic: {:x})", magic));
        };

        // Read version
        let version = read_u32(&mut reader, is_swapped)?;
        if version != 0 && version != 1 {
            return Err(anyhow!(
                "Unsupported 2bit version: {} (only version 0 and 1 are supported)",
                version
            ));
        }

        // Read seqCount
        let seq_count = read_u32(&mut reader, is_swapped)?;

        // Read reserved
        let _reserved = read_u32(&mut reader, is_swapped)?;

        // Read Index
        let mut sequence_offsets = HashMap::new();
        let mut sequence_order = Vec::new();
        for _ in 0..seq_count {
            let mut len_buf = [0u8; 1];
            reader.read_exact(&mut len_buf)?;
            let name_len = len_buf[0] as usize;

            let mut name_buf = vec![0u8; name_len];
            reader.read_exact(&mut name_buf)?;
            let name = String::from_utf8(name_buf)?;

            let offset = if version == 0 {
                read_u32(&mut reader, is_swapped)? as u64
            } else {
                read_u64(&mut reader, is_swapped)?
            };
            sequence_offsets.insert(name.clone(), offset);
            sequence_order.push(name);
        }

        Ok(Self {
            reader,
            sequence_offsets,
            sequence_order,
            is_swapped,
            version,
        })
    }

    /// Return the sequence names in their original file order.
    pub fn get_sequence_names(&self) -> Vec<String> {
        self.sequence_order.clone()
    }

    fn read_blocks(&mut self) -> Result<Blocks> {
        read_blocks(&mut self.reader, self.is_swapped)
    }

    /// Read a sequence (or subsequence) by name, applying masks unless `no_mask` is true.
    ///
    /// `start` and `end` are 0-based half-open coordinates; `None` means the full sequence.
    pub fn read_sequence(
        &mut self,
        name: &str,
        start: Option<usize>,
        end: Option<usize>,
        no_mask: bool,
    ) -> Result<String> {
        let offset = *self
            .sequence_offsets
            .get(name)
            .ok_or_else(|| anyhow!("Sequence not found: {}", name))?;

        self.reader.seek(SeekFrom::Start(offset))?;

        read_2bit_record(&mut self.reader, self.is_swapped, start, end, no_mask)
    }

    /// Return the total length (including Ns) of the named sequence.
    pub fn get_sequence_len(&mut self, name: &str) -> Result<usize> {
        let offset = *self
            .sequence_offsets
            .get(name)
            .ok_or_else(|| anyhow!("Sequence not found: {}", name))?;

        self.reader.seek(SeekFrom::Start(offset))?;
        let dna_size = read_u32(&mut self.reader, self.is_swapped)? as usize;
        Ok(dna_size)
    }

    /// Return the N-blocks and soft-mask blocks for the named sequence.
    pub fn get_sequence_blocks(&mut self, name: &str) -> Result<(Blocks, Blocks)> {
        let offset = *self
            .sequence_offsets
            .get(name)
            .ok_or_else(|| anyhow!("Sequence not found: {}", name))?;

        self.reader.seek(SeekFrom::Start(offset))?;
        let _dna_size = read_u32(&mut self.reader, self.is_swapped)?;

        let n_blocks = self.read_blocks()?;
        let mask_blocks = self.read_blocks()?;

        Ok((n_blocks, mask_blocks))
    }

    /// Return the sequence length excluding N-blocks.
    pub fn get_sequence_len_no_ns(&mut self, name: &str) -> Result<usize> {
        let offset = *self
            .sequence_offsets
            .get(name)
            .ok_or_else(|| anyhow!("Sequence not found: {}", name))?;

        self.reader.seek(SeekFrom::Start(offset))?;
        let dna_size = read_u32(&mut self.reader, self.is_swapped)? as usize;

        let n_blocks = self.read_blocks()?;

        let n_count: usize = n_blocks.iter().map(|b| b.end - b.start).sum();

        if n_count > dna_size {
            return Err(anyhow!("N blocks size {} > dna size {}", n_count, dna_size));
        }

        Ok(dna_size - n_count)
    }
}

impl<R: Read + Seek> crate::libs::io::SequenceReader for TwoBitFile<R> {
    fn read_sequence(
        &mut self,
        name: &str,
        start: Option<usize>,
        end: Option<usize>,
    ) -> anyhow::Result<String> {
        // Delegate to the inherent 5-arg form, masking kept on (no_mask=false).
        TwoBitFile::read_sequence(self, name, start, end, false)
    }
}

/// Read a little-endian u32, byte-swapping if `is_swapped` is true.
pub fn read_u32<R: Read>(reader: &mut R, is_swapped: bool) -> Result<u32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    let val = u32::from_le_bytes(buf);
    if is_swapped {
        Ok(val.swap_bytes())
    } else {
        Ok(val)
    }
}

/// Read a little-endian u64, byte-swapping if `is_swapped` is true.
pub fn read_u64<R: Read>(reader: &mut R, is_swapped: bool) -> Result<u64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    let val = u64::from_le_bytes(buf);
    if is_swapped {
        Ok(val.swap_bytes())
    } else {
        Ok(val)
    }
}

/// Read `count` little-endian u32 values, byte-swapping if `is_swapped` is true.
pub fn read_u32_vec<R: Read>(reader: &mut R, count: usize, is_swapped: bool) -> Result<Vec<u32>> {
    let mut vec = Vec::with_capacity(count);
    for _ in 0..count {
        vec.push(read_u32(reader, is_swapped)?);
    }
    Ok(vec)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn create_v1_2bit_data() -> Vec<u8> {
        let mut data = Vec::new();

        // Header (16 bytes)
        data.extend_from_slice(&TWOBIT_MAGIC.to_le_bytes()); // Magic
        data.extend_from_slice(&1u32.to_le_bytes()); // Version 1
        data.extend_from_slice(&1u32.to_le_bytes()); // SeqCount 1
        data.extend_from_slice(&0u32.to_le_bytes()); // Reserved

        // Index
        // Name: "seq1"
        let name = "seq1";
        data.push(name.len() as u8);
        data.extend_from_slice(name.as_bytes());

        // Offset calculation:
        // Header (16) + NameLen(1) + Name(4) + Offset(8) = 29 bytes
        let offset = 29u64;
        data.extend_from_slice(&offset.to_le_bytes());

        // Data Record at offset 29
        // Sequence: TCAG (4 bp)
        // T=00, C=01, A=10, G=11 -> 00011011 = 0x1B
        let dna_size = 4u32;
        data.extend_from_slice(&dna_size.to_le_bytes());

        // N Blocks
        data.extend_from_slice(&0u32.to_le_bytes()); // Count

        // Mask Blocks
        data.extend_from_slice(&0u32.to_le_bytes()); // Count

        // Reserved
        data.extend_from_slice(&0u32.to_le_bytes());

        // Packed DNA
        data.push(0x1B); // TCAG

        data
    }

    #[test]
    fn test_read_v1_basic() -> Result<()> {
        let data = create_v1_2bit_data();
        let cursor = Cursor::new(data);
        let mut tb = TwoBitFile::new(cursor)?;

        assert_eq!(tb.version, 1);

        let names = tb.get_sequence_names();
        assert_eq!(names, vec!["seq1"]);

        // Read "seq1"
        let seq = tb.read_sequence("seq1", None, None, false)?;
        assert_eq!(seq, "TCAG");

        Ok(())
    }

    fn create_v0_2bit_data() -> Vec<u8> {
        let mut data = Vec::new();

        // Header (16 bytes)
        data.extend_from_slice(&TWOBIT_MAGIC.to_le_bytes()); // Magic
        data.extend_from_slice(&0u32.to_le_bytes()); // Version 0
        data.extend_from_slice(&1u32.to_le_bytes()); // SeqCount 1
        data.extend_from_slice(&0u32.to_le_bytes()); // Reserved

        // Index
        // Name: "seq1"
        let name = "seq1";
        data.push(name.len() as u8);
        data.extend_from_slice(name.as_bytes());

        // Offset calculation:
        // Header (16) + NameLen(1) + Name(4) + Offset(4) = 25 bytes
        let offset = 25u32;
        data.extend_from_slice(&offset.to_le_bytes());

        // Data Record at offset 25
        // Sequence: TCAG (4 bp)
        // T=00, C=01, A=10, G=11 -> 00011011 = 0x1B
        let dna_size = 4u32;
        data.extend_from_slice(&dna_size.to_le_bytes());

        // N Blocks
        data.extend_from_slice(&0u32.to_le_bytes()); // Count

        // Mask Blocks
        data.extend_from_slice(&0u32.to_le_bytes()); // Count

        // Reserved
        data.extend_from_slice(&0u32.to_le_bytes());

        // Packed DNA
        data.push(0x1B); // TCAG

        data
    }

    #[test]
    fn test_read_v0_basic() -> Result<()> {
        let data = create_v0_2bit_data();
        let cursor = Cursor::new(data);
        let mut tb = TwoBitFile::new(cursor)?;

        assert_eq!(tb.version, 0);

        let names = tb.get_sequence_names();
        assert_eq!(names, vec!["seq1"]);

        // Read "seq1"
        let seq = tb.read_sequence("seq1", None, None, false)?;
        assert_eq!(seq, "TCAG");

        Ok(())
    }

    #[test]
    fn test_version_check() {
        // Construct a version 2 file (unsupported)
        let mut data = Vec::new();
        data.extend_from_slice(&TWOBIT_MAGIC.to_le_bytes());
        data.extend_from_slice(&2u32.to_le_bytes()); // Version 2
                                                     // ... rest doesn't matter as it fails fast

        let cursor = Cursor::new(data);
        let res = TwoBitFile::new(cursor);
        assert!(res.is_err());
        assert!(res
            .unwrap_err()
            .to_string()
            .contains("Unsupported 2bit version: 2"));
    }

    #[test]
    fn test_blocks_from_dna() {
        let dna = "acGTNn";
        // a=2(10), c=1(01), G=3(11), T=0(00), N=0(00), n=0(00)
        // 10 01 11 00 -> 0x9C
        // 00 00 -> 0x00

        let (packed, n_blocks, mask_blocks, size) = Blocks::from_dna(dna, true);
        assert_eq!(size, 6);
        assert_eq!(packed, vec![0x9C, 0x00]);

        // N blocks: 4..6 (N, n)
        assert_eq!(n_blocks.0, vec![4..6]);

        // Mask blocks: 0..2 (a, c), 5..6 (n)
        assert_eq!(mask_blocks.0, vec![0..2, 5..6]);
    }

    #[test]
    fn test_write_read_roundtrip() -> Result<()> {
        let seqs = vec![("seq1", "TCAG"), ("seq2", "aaNgg")];

        let mut buffer = Cursor::new(Vec::new());
        let mut writer = TwoBitWriter::new(&mut buffer);
        writer.write(&seqs, true)?;

        buffer.set_position(0);
        let mut reader = TwoBitFile::new(buffer)?;

        let mut names = reader.get_sequence_names();
        names.sort();
        assert_eq!(names, vec!["seq1", "seq2"]);

        let s1 = reader.read_sequence("seq1", None, None, false)?;
        assert_eq!(s1, "TCAG");

        let s2 = reader.read_sequence("seq2", None, None, false)?;
        assert_eq!(s2, "aaNgg");

        Ok(())
    }

    // --- Independent tests for shared read_2bit_record / write_2bit_record ---

    #[test]
    fn test_record_roundtrip_acgt() -> Result<()> {
        let dna = "TCAGTCAG";
        let mut buf = Vec::new();
        write_2bit_record(&mut buf, dna, false)?;

        let mut cursor = Cursor::new(buf);
        let seq = read_2bit_record(&mut cursor, false, None, None, false)?;
        assert_eq!(seq, dna);
        Ok(())
    }

    #[test]
    fn test_record_roundtrip_with_mask() -> Result<()> {
        let dna = "aaGTcc"; // lowercase = soft mask
        let mut buf = Vec::new();
        write_2bit_record(&mut buf, dna, true)?;

        let mut cursor = Cursor::new(buf.clone());
        // with mask applied
        let seq = read_2bit_record(&mut cursor, false, None, None, false)?;
        assert_eq!(seq, dna);

        let mut cursor = Cursor::new(buf);
        // no mask -> uppercase
        let seq_no_mask = read_2bit_record(&mut cursor, false, None, None, true)?;
        assert_eq!(seq_no_mask, "AAGTCC");
        Ok(())
    }

    #[test]
    fn test_record_roundtrip_with_n() -> Result<()> {
        let dna = "ACGTNNNNTTTT";
        let mut buf = Vec::new();
        write_2bit_record(&mut buf, dna, false)?;

        let mut cursor = Cursor::new(buf);
        let seq = read_2bit_record(&mut cursor, false, None, None, false)?;
        assert_eq!(seq, dna);
        Ok(())
    }

    #[test]
    fn test_record_slice() -> Result<()> {
        let dna = "ACGTACGTACGT";
        let mut buf = Vec::new();
        write_2bit_record(&mut buf, dna, false)?;

        let mut cursor = Cursor::new(buf);
        let seq = read_2bit_record(&mut cursor, false, Some(2), Some(8), false)?;
        assert_eq!(seq, &dna[2..8]);
        Ok(())
    }

    #[test]
    fn test_record_empty_slice() -> Result<()> {
        let dna = "ACGT";
        let mut buf = Vec::new();
        write_2bit_record(&mut buf, dna, false)?;

        let mut cursor = Cursor::new(buf);
        let seq = read_2bit_record(&mut cursor, false, Some(2), Some(2), false)?;
        assert_eq!(seq, "");
        Ok(())
    }

    #[test]
    fn test_record_non_acgt_becomes_n() -> Result<()> {
        // U/R and IUPAC ambiguity codes are not valid DNA bases; they should be
        // stored as N-blocks rather than silently packed as T.
        let dna = "ACGTURacgtur";
        let mut buf = Vec::new();
        write_2bit_record(&mut buf, dna, true)?;

        let mut cursor = Cursor::new(buf.clone());
        let seq = read_2bit_record(&mut cursor, false, None, None, false)?;
        assert_eq!(seq, "ACGTNNacgtnn");

        let mut cursor = Cursor::new(buf);
        let seq_no_mask = read_2bit_record(&mut cursor, false, None, None, true)?;
        assert_eq!(seq_no_mask, "ACGTNNACGTNN");

        Ok(())
    }
}
