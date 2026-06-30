//! FASTA helpers wrapping noodles_fasta with pgr I/O.
//!
//! Division of labor for FASTA-related code:
//! * [`fmt::fa`] (this module) — sequential read/write of FASTA streams,
//!   record construction, windowing.
//! * [`libs::loc`] — random-access sequence extraction by genomic interval
//!   (uses 2bit/BGZF FastaStore backends).
//! * [`libs::fasta::stat`] — sequence statistics (N50, base counts, etc.).

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

/// Run the `fa window` workflow: split sequences into overlapping windows
/// and optionally chunk/shuffle the output records.
///
/// * `infile` — input FASTA path (supports stdin/.gz via [`io::reader`]).
/// * `len` / `step` — window length and step size in bases.
/// * `shuffle` — randomize record order (uses `seed` for reproducibility).
/// * `chunk_size` — when set, split output into files of N records each
///   (`outfile` must not be `stdout`).
/// * `outfile` — `stdout` or a file path; chunked files are named
///   `<stem>.NNN<ext>`.
///
/// Windows consisting entirely of N/n are skipped (see [`windows`]).
pub fn run_window(
    infile: &str,
    len: usize,
    step: usize,
    shuffle: bool,
    seed: u64,
    chunk_size: Option<usize>,
    outfile: &str,
) -> anyhow::Result<()> {
    if chunk_size.is_some() && outfile == "stdout" {
        anyhow::bail!("Cannot use --chunk with stdout output");
    }

    let mut fa_in = reader(infile)?;

    // Build a chunked output path: <stem>.NNN<ext>
    let create_writer = |part: usize| -> anyhow::Result<Box<dyn std::io::Write>> {
        if outfile == "stdout" {
            io::writer("stdout")
        } else {
            let path = std::path::Path::new(outfile);
            let file_stem = path
                .file_stem()
                .and_then(std::ffi::OsStr::to_str)
                .ok_or_else(|| anyhow::anyhow!("invalid outfile stem: {}", outfile))?;
            let extension = path
                .extension()
                .and_then(std::ffi::OsStr::to_str)
                .unwrap_or_default();
            let ext_str = if extension.is_empty() {
                String::new()
            } else {
                format!(".{}", extension)
            };
            let new_filename = format!("{}.{:03}{}", file_stem, part, ext_str);
            let new_path = path.with_file_name(new_filename);
            let new_path_str = new_path
                .to_str()
                .ok_or_else(|| anyhow::anyhow!("invalid chunked path: {}", new_path.display()))?;
            io::writer(new_path_str)
        }
    };

    let mut current_part = 1;
    let mut record_count = 0;

    let mut fa_out: Option<fasta::io::Writer<Box<dyn std::io::Write>>> = None;

    // Initialize global writer if not chunking.
    if chunk_size.is_none() {
        fa_out = Some(writer(outfile)?);
    } else if !shuffle {
        // If chunking without shuffle, init first writer
        let w = create_writer(current_part)?;
        fa_out = Some(writer_from_writer(w));
    }

    // For shuffle we accumulate records; for non-shuffle chunking we stream.
    let mut records_buffer: Vec<fasta::Record> = Vec::new();

    for result in fa_in.records() {
        let record = result?;
        let name = String::from_utf8(record.name().into())?;
        let seq = record.sequence();

        for (new_name, window) in windows(&name, seq.as_ref(), len, step) {
            let new_record = new_record(&new_name, &window);

            if shuffle {
                records_buffer.push(new_record);

                // If chunk limit reached, flush buffer
                if let Some(limit) = chunk_size {
                    if records_buffer.len() >= limit {
                        flush_shuffled_chunk(
                            &mut records_buffer,
                            seed,
                            current_part,
                            &create_writer,
                        )?;
                        current_part += 1;
                    }
                }
            } else {
                // No shuffle
                if let Some(limit) = chunk_size {
                    if record_count >= limit {
                        current_part += 1;
                        record_count = 0;
                        let w = create_writer(current_part)?;
                        fa_out = Some(writer_from_writer(w));
                    }
                }

                if let Some(ref mut w) = fa_out {
                    w.write_record(&new_record)?;
                    record_count += 1;
                }
            }
        }
    }

    // Flush remaining buffer (Shuffle case)
    if shuffle && !records_buffer.is_empty() {
        use rand::seq::SliceRandom;
        use rand::SeedableRng;
        let chunk_seed = seed + (current_part as u64);
        let mut rng = rand::rngs::StdRng::seed_from_u64(chunk_seed);
        records_buffer.shuffle(&mut rng);

        // If chunking, this goes to a new chunk file.
        // If not chunking, this goes to the single global file.
        let mut final_out = if chunk_size.is_some() {
            let w = create_writer(current_part)?;
            writer_from_writer(w)
        } else if let Some(w) = fa_out.take() {
            w
        } else {
            writer(outfile)?
        };

        for record in records_buffer {
            final_out.write_record(&record)?;
        }
    }

    Ok(())
}

// Helper: shuffle `records_buffer` with a per-chunk seed and write to the
// chunk file identified by `part`. Clears the buffer on success.
fn flush_shuffled_chunk(
    records_buffer: &mut Vec<fasta::Record>,
    seed: u64,
    part: usize,
    create_writer: &impl Fn(usize) -> anyhow::Result<Box<dyn std::io::Write>>,
) -> anyhow::Result<()> {
    use rand::seq::SliceRandom;
    use rand::SeedableRng;
    let chunk_seed = seed + (part as u64);
    let mut rng = rand::rngs::StdRng::seed_from_u64(chunk_seed);
    records_buffer.shuffle(&mut rng);

    let w = create_writer(part)?;
    let mut chunk_out = writer_from_writer(w);
    for r in records_buffer.iter() {
        chunk_out.write_record(r)?;
    }
    records_buffer.clear();
    Ok(())
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
