//! FASTA helpers: owned records and a self-contained writer (no noodles).
//!
//! Division of labor for FASTA-related code:
//! * [`fmt::fa`] (this module) — sequential read/write of FASTA streams,
//!   record construction, windowing.
//! * [`libs::loc`] — random-access sequence extraction by genomic interval
//!   (uses 2bit/BGZF FastaStore backends).
//! * [`libs::fasta::stat`] — sequence statistics (N50, base counts, etc.).

use crate::libs::fmt::seq::{SeqReader, SeqRecord};
use crate::libs::io as pgr_io;
use bstr::{BStr, BString};
use std::io::{self, Write};

/// Owned FASTA record: name, optional description, sequence.
#[derive(Debug, Clone, Default)]
pub struct FastaRecord {
    name: BString,
    desc: Option<BString>,
    seq: Vec<u8>,
}

impl FastaRecord {
    /// Builds a record from a name and a sequence byte slice.
    pub fn new(name: &str, seq: &[u8]) -> Self {
        Self {
            name: BString::from(name),
            desc: None,
            seq: seq.to_vec(),
        }
    }

    /// Builds a record from a name, optional description, and sequence.
    pub fn with_desc(name: &str, desc: Option<&[u8]>, seq: &[u8]) -> Self {
        Self {
            name: BString::from(name),
            desc: desc.map(BString::from),
            seq: seq.to_vec(),
        }
    }

    /// Builds a record from a new name and sequence, preserving the
    /// description from `source`.
    pub fn preserving_desc(name: &str, source: &Self, seq: &[u8]) -> Self {
        Self {
            name: BString::from(name),
            desc: source.desc.clone(),
            seq: seq.to_vec(),
        }
    }

    /// Record name.
    pub fn name(&self) -> &BStr {
        self.name.as_ref()
    }

    /// Optional description (text after the name).
    pub fn description(&self) -> Option<&[u8]> {
        self.desc.as_deref().map(|d| d.as_ref())
    }

    /// Sequence bases.
    pub fn sequence(&self) -> &[u8] {
        &self.seq
    }
}

/// A self-contained FASTA writer.
///
/// `line_base_count` is the number of bases per line; `usize::MAX` writes
/// each sequence on a single line (the default).
pub struct FastaWriter<W: Write> {
    inner: W,
    line_base_count: usize,
}

impl<W: Write> FastaWriter<W> {
    /// Creates a writer over `inner` with no line wrapping.
    pub fn new(inner: W) -> Self {
        Self {
            inner,
            line_base_count: usize::MAX,
        }
    }

    /// Creates a writer with a configurable line width.
    pub fn with_line_width(inner: W, line_base_count: usize) -> Self {
        Self {
            inner,
            line_base_count,
        }
    }

    /// Writes one record (`>name[ desc]\nseq\n`).
    pub fn write_record(&mut self, record: &FastaRecord) -> io::Result<()> {
        write!(self.inner, ">{}", record.name)?;
        if let Some(desc) = &record.desc {
            write!(self.inner, " {desc}")?;
        }
        self.inner.write_all(b"\n")?;
        if !record.seq.is_empty() {
            write_seq(&mut self.inner, &record.seq, self.line_base_count)?;
        }
        Ok(())
    }

    /// Access to the underlying writer.
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Consumes the writer, returning the underlying writer.
    pub fn into_inner(self) -> W {
        self.inner
    }
}

fn write_seq<W: Write>(w: &mut W, seq: &[u8], line_base_count: usize) -> io::Result<()> {
    if seq.len() <= line_base_count {
        w.write_all(seq)?;
        w.write_all(b"\n")?;
    } else {
        for chunk in seq.chunks(line_base_count) {
            w.write_all(chunk)?;
            w.write_all(b"\n")?;
        }
    }
    Ok(())
}

/// Open a FASTA reader from a path (supports stdin and gzip).
pub fn reader(infile: &str) -> anyhow::Result<SeqReader<'static>> {
    SeqReader::new(infile)
}

/// Create a FASTA writer with no line wrapping (single-line sequences).
pub fn writer(outfile: &str) -> anyhow::Result<FastaWriter<Box<dyn std::io::Write>>> {
    Ok(writer_from_writer(Box::new(pgr_io::writer(outfile)?)))
}

/// Create a FASTA writer with configurable line width.
pub fn writer_with_wrap(
    outfile: &str,
    line_base_count: usize,
) -> anyhow::Result<FastaWriter<Box<dyn std::io::Write>>> {
    let w: Box<dyn std::io::Write> = Box::new(pgr_io::writer(outfile)?);
    Ok(FastaWriter::with_line_width(w, line_base_count))
}

/// Wrap an existing writer as a FASTA writer with no line wrapping.
pub fn writer_from_writer<W: std::io::Write>(writer: W) -> FastaWriter<W> {
    FastaWriter::new(writer)
}

/// Build a FASTA record from a name and a sequence byte slice.
pub fn new_record(name: &str, seq: &[u8]) -> FastaRecord {
    FastaRecord::new(name, seq)
}

/// Build a FASTA record from a name, optional description, and sequence.
pub fn new_record_with_desc(name: &str, desc: Option<&[u8]>, seq: &[u8]) -> FastaRecord {
    FastaRecord::with_desc(name, desc, seq)
}

/// Build a FASTA record from a new name and sequence, preserving the
/// description from `source`.
pub fn new_record_preserving_desc(name: &str, source: &FastaRecord, seq: &[u8]) -> FastaRecord {
    FastaRecord::preserving_desc(name, source, seq)
}

/// Generate windowed sub-sequences from `name`/`seq`.
///
/// Each window is `len` bytes long (the last one may be shorter); successive
/// windows start `step` bytes apart. Windows consisting entirely of ambiguous
/// bases (N/IUPAC) are skipped. Coordinates embedded in the emitted names are
/// 1-based inclusive (`name:start-end`).
pub fn windows(name: &str, seq: &[u8], len: usize, step: usize) -> Vec<(String, Vec<u8>)> {
    if step == 0 {
        return Vec::new();
    }
    let mut result = Vec::new();
    let seq_len = seq.len();
    for start in (0..seq_len).step_by(step) {
        // `saturating_add` guards against a window length large enough that
        // `start + len` would overflow `usize` (e.g. `--window` near the max).
        let end = std::cmp::min(start.saturating_add(len), seq_len);
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
/// Windows consisting entirely of ambiguous bases (N/IUPAC) are skipped (see
/// [`windows`]).
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
        anyhow::bail!("Cannot use --chunk-records with stdout output");
    }

    let mut fa_in = reader(infile)?;
    let mut rec = SeqRecord::default();

    // Build a chunked output path: <stem>.NNN<ext>
    let create_writer = |part: usize| -> anyhow::Result<Box<dyn std::io::Write>> {
        if outfile == "stdout" {
            Ok(Box::new(pgr_io::writer("stdout")?))
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
            // A chunk filename derived from `-o` (e.g. `out.001.fa`) can
            // collide with the input file when the input happens to share that
            // name; opening it with `truncate` would clobber the input while
            // it is still being streamed. Reject before touching the file.
            if chunk_size.is_some() && crate::libs::io::same_path(new_path_str, infile) {
                anyhow::bail!(
                    "chunked output file {} would overwrite input file {}",
                    new_path_str,
                    infile
                );
            }
            Ok(Box::new(pgr_io::writer(new_path_str)?))
        }
    };

    let mut current_part = 1;
    let mut record_count = 0;

    let mut fa_out: Option<FastaWriter<Box<dyn std::io::Write>>> = None;

    // Initialize global writer if not chunking.
    if chunk_size.is_none() {
        fa_out = Some(writer(outfile)?);
    } else if !shuffle {
        // If chunking without shuffle, init first writer
        let w = create_writer(current_part)?;
        fa_out = Some(writer_from_writer(w));
    }

    // For shuffle we accumulate records; for non-shuffle chunking we stream.
    let mut records_buffer: Vec<FastaRecord> = Vec::new();

    while fa_in.read_record(&mut rec)? {
        let name = String::from_utf8(rec.name().to_vec())?;
        let seq = rec.sequence();

        for (new_name, window) in windows(&name, seq, len, step) {
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
                        if let Some(ref mut w) = fa_out {
                            w.get_mut().flush()?;
                        }
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
        final_out.get_mut().flush()?;
    }

    // Flush streaming writer (non-shuffle path)
    if let Some(ref mut w) = fa_out {
        w.get_mut().flush()?;
    }

    Ok(())
}

// Helper: shuffle `records_buffer` with a per-chunk seed and write to the
// chunk file identified by `part`. Clears the buffer on success.
fn flush_shuffled_chunk(
    records_buffer: &mut Vec<FastaRecord>,
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
    chunk_out.get_mut().flush()?;
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
    crate::libs::bgzf::build_gzi_index(path)
}

/// Recursively collect FASTA files (`.fa` and `.fa.gz`) under `path`.
/// A file input is returned as a single-element vec. Directory inputs are
/// walked recursively, matching `.fa` and `.fa.gz` extensions.
pub fn find_fasta_files<P: AsRef<std::path::Path>>(path: P) -> Vec<std::path::PathBuf> {
    let mut files = Vec::new();
    let path = path.as_ref();

    if path.is_file() {
        files.push(path.to_path_buf());
    } else if path.is_dir() {
        if let Ok(entries) = std::fs::read_dir(path) {
            for entry in entries.flatten() {
                let p = entry.path();
                if p.is_dir() {
                    files.extend(find_fasta_files(&p));
                } else if let Some(ext) = p.extension() {
                    let ext_str = ext.to_string_lossy().to_lowercase();
                    if ext_str == "fa" {
                        files.push(p);
                    } else if ext_str == "gz" {
                        if let Some(stem) = p.file_stem() {
                            let stem_path = std::path::Path::new(stem);
                            if let Some(stem_ext) = stem_path.extension() {
                                if stem_ext.to_string_lossy().to_lowercase() == "fa" {
                                    files.push(p);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    files
}

/// Mask sequence regions. Soft-mask lowercases, hard-mask replaces with N.
///
/// Operates on raw bytes so mask spans are interpreted as byte offsets; this
/// also avoids panicking on char boundaries when a sequence contains
/// multi-byte UTF-8 (a `&str` slice would panic mid-character).
pub fn mask_sequence(
    seq: &[u8],
    spans: &crate::libs::ds::IntSpan,
    hard: bool,
) -> anyhow::Result<Vec<u8>> {
    let mut out = seq.to_vec();
    for (lower, upper) in spans.spans().iter() {
        if *lower < 1 {
            anyhow::bail!("span lower bound must be >= 1, got {}", lower);
        }
        let offset = (*lower - 1) as usize;
        let length = (*upper - *lower + 1) as usize;
        if offset + length > out.len() {
            anyhow::bail!(
                "span {}-{} exceeds sequence length {}",
                lower,
                upper,
                out.len()
            );
        }
        if hard {
            out[offset..offset + length].fill(b'N');
        } else {
            out[offset..offset + length].make_ascii_lowercase();
        }
    }
    Ok(out)
}

/// Find contiguous masked regions (lowercase and/or N/n) in a sequence. Returns 0-based inclusive (begin, end) pairs.
pub fn find_masked_regions(seq: &[u8], gap_only: bool) -> Vec<(usize, usize)> {
    let words = crate::libs::nt_simd::masked_bitmap(seq, gap_only);
    let mut regions = Vec::new();
    let mut begin: Option<usize> = None;
    let mut end: Option<usize> = None;

    for (wi, word) in words.into_iter().enumerate() {
        let base = wi * 32;
        let mut w = word;
        let mut off = 0usize;
        while w != 0 {
            let tz = w.trailing_zeros() as usize;
            let run = w >> tz;
            let ones = run.trailing_ones() as usize;
            let s = base + off + tz;
            let e = s + ones - 1;
            match end {
                Some(prev_end) if s == prev_end + 1 => end = Some(e),
                Some(prev_end) => {
                    regions.push((begin.unwrap(), prev_end));
                    begin = Some(s);
                    end = Some(e);
                }
                None => {
                    begin = Some(s);
                    end = Some(e);
                }
            }
            off += tz + ones;
            w >>= tz + ones;
        }
    }
    if let (Some(b), Some(e)) = (begin, end) {
        regions.push((b, e));
    }
    regions
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn windows_huge_length_does_not_overflow() {
        // Regression: `start + len` used to overflow `usize` when `len` was
        // near the max and a multi-window sequence produced a large `start`.
        let seq = vec![b'A'; 200];
        let windows = windows("seq", &seq, usize::MAX - 100, 1);
        // One window per start position (window length far exceeds the
        // sequence, so the window is always the whole remaining suffix).
        assert_eq!(windows.len(), 200);
    }

    #[test]
    fn windows_empty_and_zero_step() {
        assert!(windows("seq", b"", 10, 5).is_empty());
        assert!(windows("seq", b"ACGT", 10, 0).is_empty());
    }

    #[test]
    fn mask_sequence_non_ascii_does_not_panic() {
        // Regression: a multi-byte UTF-8 base used to panic on a `&str` byte
        // slice when the mask span landed mid-character.
        let seq = b"A\xc3\xa9C"; // A, é (2 bytes), C
        let spans = crate::libs::ds::IntSpan::from_pair(2, 2);
        // Soft mask: byte 2 = 0xc3; ASCII lowercase leaves it unchanged.
        let masked = mask_sequence(seq, &spans, false).unwrap();
        assert_eq!(masked, seq.to_vec());
        // Hard mask: byte 2 (0xc3) is replaced with 'N'.
        let hard = mask_sequence(seq, &spans, true).unwrap();
        assert_eq!(hard, b"A\x4e\xa9C".to_vec());
    }
}
