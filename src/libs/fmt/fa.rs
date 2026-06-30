//! FASTA helpers wrapping noodles_fasta with pgr I/O.

use crate::libs::io;
use noodles_bgzf as bgzf;
use noodles_fasta as fasta;
use std::io::{Read, Seek};

/// Open a FASTA reader from a path (supports stdin and gzip).
pub fn reader(infile: &str) -> anyhow::Result<fasta::io::Reader<Box<dyn std::io::BufRead>>> {
    let r = io::reader(infile)?;
    Ok(fasta::io::Reader::new(r))
}

/// Create a FASTA writer with no line wrapping (single-line sequences).
pub fn writer(outfile: &str) -> anyhow::Result<fasta::io::Writer<Box<dyn std::io::Write>>> {
    Ok(writer_from_writer(io::writer(outfile)?))
}

/// Create a FASTA writer with configurable line width.
pub fn writer_with_wrap(
    outfile: &str,
    line_base_count: usize,
) -> anyhow::Result<fasta::io::Writer<Box<dyn std::io::Write>>> {
    let w = io::writer(outfile)?;
    Ok(fasta::io::writer::Builder::default()
        .set_line_base_count(line_base_count)
        .build_from_writer(w))
}

/// Wrap an existing writer as a FASTA writer with no line wrapping.
pub fn writer_from_writer<W: std::io::Write>(writer: W) -> fasta::io::Writer<W> {
    fasta::io::writer::Builder::default()
        .set_line_base_count(usize::MAX)
        .build_from_writer(writer)
}

/// Build a FASTA record from a name and a sequence byte slice.
pub fn new_record(name: &str, seq: &[u8]) -> fasta::Record {
    let definition = fasta::record::Definition::new(name, None);
    let sequence = fasta::record::Sequence::from(seq.to_vec());
    fasta::Record::new(definition, sequence)
}

/// Generate windowed sub-sequences from `name`/`seq`.
///
/// Each window is `len` bytes long (the last one may be shorter); successive
/// windows start `step` bytes apart. Windows consisting entirely of N/n are
/// skipped. Coordinates embedded in the emitted names are 1-based inclusive
/// (`name:start-end`).
pub fn windows(name: &str, seq: &[u8], len: usize, step: usize) -> Vec<(String, Vec<u8>)> {
    let mut result = Vec::new();
    let seq_len = seq.len();
    for start in (0..seq_len).step_by(step) {
        let end = std::cmp::min(start + len, seq_len);
        if start >= end {
            continue;
        }
        let window = &seq[start..end];
        if window.iter().all(|&b| crate::libs::nt::is_n(b)) {
            continue;
        }
        let new_name = format!("{}:{}-{}", name, start + 1, end);
        result.push((new_name, window.to_vec()));
    }
    result
}

/// Build a .gzi index for a BGZF file.
///
/// The GZI format is defined by `bgzip` and used for random access.
/// It consists of:
/// 1. A header (u64): Number of entries
/// 2. Entries (pairs of u64): (compressed_offset, uncompressed_offset)
///
/// Note:
/// * The format is Little-Endian.
/// * The first BGZF block (offset 0, 0) is implicitly skipped and NOT included in the index.
/// * Empty blocks (like EOF markers with ISIZE=0) are also skipped.
pub fn build_gzi_index(path: &str) -> anyhow::Result<()> {
    let mut file = std::fs::File::open(path)?;
    let mut index_data = Vec::new();
    let mut uncompressed_offset = 0;
    let mut compressed_offset = 0;

    loop {
        file.seek(std::io::SeekFrom::Start(compressed_offset))?;

        let mut header_fixed = [0u8; 12];
        match file.read_exact(&mut header_fixed) {
            Ok(_) => {}
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
            Err(e) => return Err(e.into()),
        }

        if header_fixed[0] != 0x1f || header_fixed[1] != 0x8b {
            break;
        }

        let flg = header_fixed[3];
        if (flg & 4) == 0 {
            break;
        }

        let xlen = u16::from_le_bytes([header_fixed[10], header_fixed[11]]) as u64;
        if xlen == 0 {
            break;
        }

        let mut extra = vec![0u8; xlen as usize];
        file.read_exact(&mut extra)?;

        let mut bsize = 0u16;
        let mut cursor = 0;
        let mut found_bc = false;

        while cursor + 4 <= extra.len() {
            let si1 = extra[cursor];
            let si2 = extra[cursor + 1];
            let slen = u16::from_le_bytes([extra[cursor + 2], extra[cursor + 3]]);

            if si1 == b'B' && si2 == b'C' && slen == 2 {
                if cursor + 6 <= extra.len() {
                    bsize = u16::from_le_bytes([extra[cursor + 4], extra[cursor + 5]]);
                    found_bc = true;
                }
                break;
            }
            cursor += 4 + slen as usize;
        }

        if !found_bc {
            return Err(anyhow::anyhow!(
                "Missing BC subfield in BGZF block at offset {}",
                compressed_offset
            ));
        }

        let block_size = bsize as u64 + 1;

        file.seek(std::io::SeekFrom::Start(compressed_offset + block_size - 4))?;
        let mut isize_buf = [0u8; 4];
        file.read_exact(&mut isize_buf)?;
        let isize = u32::from_le_bytes(isize_buf) as u64;

        if compressed_offset > 0 && isize > 0 {
            index_data.push((compressed_offset, uncompressed_offset));
        }

        compressed_offset += block_size;
        uncompressed_offset += isize;
    }

    let index = bgzf::gzi::Index::from(index_data);
    let index_path = format!("{}.gzi", path);
    let mut writer = std::fs::File::create(index_path)?;
    bgzf::gzi::io::Writer::new(&mut writer).write_index(&index)?;

    Ok(())
}
