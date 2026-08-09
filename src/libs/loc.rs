use crate::libs::bgzf::CachedBgzfReader;
use crate::libs::ds::Range;
use crate::libs::fmt::fa::FastaRecord;
use crate::libs::fmt::seq::{SeqReader, SeqRecord};
use indexmap::IndexMap;
use std::io::BufReader;
use std::io::{Read, Seek, SeekFrom, Write};
use std::num::NonZeroUsize;

/// Random-access reader for indexed FASTA files (plain or BGZF-compressed).
pub enum Input {
    File(std::fs::File),
    Bgzf(Box<CachedBgzfReader>),
}

/// Default number of decompressed BGZF blocks (<= 64 KiB each) to cache.
const BGZF_BLOCK_CACHE: usize = 16;

pub fn create_loc(infile: &str, locfile: &str, is_bgzf: bool) -> anyhow::Result<()> {
    let mut reader: Box<dyn std::io::BufRead> = if is_bgzf {
        // http://www.htslib.org/doc/bgzip.html
        // Bgzip will attempt to ensure BGZF blocks end on a newline when the input is a text file.
        // The exception to this is where a single line is larger than a BGZF block (64Kb).
        Box::new(BufReader::new(crate::libs::bgzf::GzReader::new(
            std::fs::File::open(infile)?,
        )?))
    } else {
        crate::libs::io::reader(infile)?
    };

    let mut writer = crate::libs::io::writer(locfile)?;

    // https://www.ginkgobioworks.com/2023/03/17/even-more-rapid-retrieval-from-very-large-files-with-rust/
    let mut record_size = 0; // including header, sequence, newlines
    let mut offset = 0;
    let mut line = String::new();
    loop {
        let num = reader.read_line(&mut line)?;
        if num == 0 {
            break;
        }

        if let Some(stripped) = line.strip_prefix('>') {
            if record_size > 0 {
                // the size of the previous record
                writer.write_fmt(format_args!("\t{}\n", record_size))?;
            }
            // reset size counter for new record
            record_size = 0;

            //current record name
            let name = stripped
                .split(|c: char| c.is_ascii_whitespace())
                .next()
                .unwrap_or("");
            writer.write_fmt(format_args!("{}\t{}", name, offset))?;
        }

        record_size += num;
        offset += num;
        line.clear();
    }
    if record_size > 0 {
        writer.write_fmt(format_args!("\t{}\n", record_size))?;
    }

    Ok(())
}

/// Open a FASTA file as `Input` (plain `File` or BGZF `IndexedReader`).
pub fn open_input(infile: &str, is_bgzf: bool) -> anyhow::Result<Input> {
    if is_bgzf {
        let capacity = NonZeroUsize::new(BGZF_BLOCK_CACHE).expect("non-zero cache size");
        Ok(Input::Bgzf(Box::new(CachedBgzfReader::open(
            infile, capacity,
        )?)))
    } else {
        Ok(Input::File(std::fs::File::open(std::path::Path::new(
            infile,
        ))?))
    }
}

/// Open a FASTA file with .loc index for random access.
/// Creates the .loc index if it doesn't exist, if `force_update` is true, or
/// if the existing index is older than the FASTA file (stale index).
/// Returns the Input reader and the loaded .loc index.
#[allow(clippy::type_complexity)]
pub fn open_indexed(
    infile: &str,
    force_update: bool,
) -> anyhow::Result<(Input, IndexMap<String, (u64, usize)>)> {
    let is_bgzf = crate::is_bgzf(infile);
    let loc_file = format!("{}.loc", infile);
    if !std::path::Path::new(&loc_file).is_file()
        || force_update
        || !loc_is_fresh(infile, &loc_file)
    {
        create_loc(infile, &loc_file, is_bgzf)?;
    }
    let loc_of = load_loc(&loc_file)?;
    let reader = open_input(infile, is_bgzf)?;
    Ok((reader, loc_of))
}

/// True when the `.loc` index exists and is not older than the FASTA file.
///
/// A stale index (e.g. the FASTA was edited after the index was built) would
/// serve wrong offsets/sizes, so callers rebuild it instead.
fn loc_is_fresh(infile: &str, loc_file: &str) -> bool {
    let Ok(fa_meta) = std::fs::metadata(infile) else {
        return false;
    };
    let Ok(loc_meta) = std::fs::metadata(loc_file) else {
        return false;
    };
    if !loc_meta.is_file() {
        return false;
    }
    match (fa_meta.modified(), loc_meta.modified()) {
        (Ok(fa_m), Ok(loc_m)) => loc_m >= fa_m,
        // mtimes unavailable (e.g. unusual filesystems): keep the existing
        // index rather than rebuilding on every call.
        _ => true,
    }
}

pub fn load_loc(loc_file: &str) -> anyhow::Result<IndexMap<String, (u64, usize)>> {
    let mut reader = crate::libs::io::reader(loc_file)?;

    let mut loc_of: IndexMap<String, (u64, usize)> = IndexMap::new();
    let mut line = String::new();
    while let Ok(num) = reader.by_ref().read_line(&mut line) {
        if num == 0 {
            break;
        }
        let fields: Vec<&str> = line.trim().split('\t').collect();
        if fields.len() != 3 {
            continue;
        }

        loc_of.insert(
            fields[0].to_string(),
            (fields[1].parse::<u64>()?, fields[2].parse::<usize>()?),
        );

        line.clear();
    }

    Ok(loc_of)
}

pub fn fetch_record(
    reader: &mut Input,
    loc_of: &IndexMap<String, (u64, usize)>,
    name: &str,
) -> anyhow::Result<FastaRecord> {
    let (offset, size) = loc_of
        .get(name)
        .ok_or_else(|| anyhow::anyhow!("{} not found in the .loc index file", name))?;

    let data_buf = read_offset(reader, *offset, *size)?;
    let mut fa_in = SeqReader::from_reader(Box::new(std::io::Cursor::new(data_buf)));
    let mut rec = SeqRecord::default();
    if !fa_in.read_record(&mut rec)? {
        anyhow::bail!("empty record for {}", name);
    }
    // The .loc index stores the bare name; the description is not part of the
    // indexed record (historical behavior).
    Ok(FastaRecord::new(name, rec.sequence()))
}

pub fn records_offset(
    reader: &mut Input,
    offset: u64,
    size: usize,
) -> anyhow::Result<Vec<FastaRecord>> {
    let mut records = Vec::new();

    let data_buf = read_offset(reader, offset, size)?;
    let mut fa_in = SeqReader::from_reader(Box::new(std::io::Cursor::new(data_buf)));
    let mut rec = SeqRecord::default();
    while fa_in.read_record(&mut rec)? {
        let name = String::from_utf8_lossy(rec.name());
        records.push(FastaRecord::with_desc(
            &name,
            rec.description(),
            rec.sequence(),
        ));
    }

    Ok(records)
}

/// Slice a subsequence from `record` by 1-based `rg`, applying reverse
/// complement for `-` strand. Returns the resulting owned sequence.
pub fn slice_record(record: &FastaRecord, rg: &crate::libs::ds::Range) -> anyhow::Result<Vec<u8>> {
    let seq = record.sequence();
    let start = *rg.start() as usize;
    let end = *rg.end() as usize;
    if start == 0 || end < start || end > seq.len() {
        anyhow::bail!("slice error for [{}]", rg);
    }
    let mut slice = seq[start - 1..end].to_vec();
    if rg.strand() == "-" {
        // Reverse complement using the `NT_COMP` lookup table (standard and
        // IUPAC bases complemented, case preserved; unknown bytes like `-`/`*`
        // kept as-is). This matches `fa rc`'s documented behavior and avoids
        // `Sequence::complement()`, which errors on non-IUPAC characters.
        slice = slice
            .iter()
            .rev()
            .map(|&b| {
                let c = crate::libs::nt::NT_COMP[b as usize];
                if c == 255 {
                    b
                } else {
                    c
                }
            })
            .collect();
    }
    Ok(slice)
}

pub fn fetch_range_seq(
    reader: &mut Input,
    loc_of: &IndexMap<String, (u64, usize)>,
    rg: &crate::libs::ds::Range,
) -> anyhow::Result<String> {
    let seq_id = rg.chr();
    if !loc_of.contains_key(seq_id) {
        return Err(anyhow::anyhow!(
            "{} for [{}] not found in the .loc index file",
            seq_id,
            rg
        ));
    }

    let record = fetch_record(reader, loc_of, seq_id)?;

    // name only
    if *rg.start() == 0 {
        let seq = record.sequence().iter().map(|&b| b as char).collect();
        return Ok(seq);
    }

    let slice = slice_record(&record, rg)?;
    let seq = slice.iter().map(|&b| b as char).collect();
    Ok(seq)
}

pub fn read_offset(reader: &mut Input, offset: u64, size: usize) -> anyhow::Result<Vec<u8>> {
    let mut data_buf = vec![0; size];

    match reader {
        Input::File(rdr) => {
            rdr.seek(SeekFrom::Start(offset))?;
            rdr.read_exact(&mut data_buf)?;
        }
        Input::Bgzf(rdr) => {
            rdr.seek(SeekFrom::Start(offset))?;
            rdr.read_exact(&mut data_buf)?;
        }
    }

    Ok(data_buf)
}

/// ```
/// let seq = pgr::libs::loc::get_seq_loc("tests/fas/NC_000932.fa", "NC_000932:1-10").unwrap();
/// assert_eq!(seq, "ATGGGCGAAC".to_string());
/// let seq = pgr::libs::loc::get_seq_loc("tests/fas/NC_000932.fa", "NC_000932(-):1-10").unwrap();
/// assert_eq!(seq, "GTTCGCCCAT".to_string());
/// let res = pgr::libs::loc::get_seq_loc("tests/fas/NC_000932.fa", "FAKE:1-10");
/// assert_eq!(res.unwrap(), "".to_string());
/// ```
pub fn get_seq_loc(file: &str, range: &str) -> anyhow::Result<String> {
    let range = Range::from_str(range);
    if !range.is_valid() {
        return Ok("".to_string());
    }

    let (mut reader, loc_of) = open_indexed(file, false)?;

    if !loc_of.contains_key(range.chr()) {
        return Ok("".to_string());
    }

    let seq = fetch_range_seq(&mut reader, &loc_of, &range)?;

    Ok(seq)
}

/// Merge overlapping or adjacent half-open intervals.
pub fn merge_intervals(mut blocks: Vec<std::ops::Range<usize>>) -> Vec<std::ops::Range<usize>> {
    blocks.sort_by_key(|r| r.start);
    let mut merged: Vec<std::ops::Range<usize>> = Vec::new();
    if let Some(first) = blocks.first() {
        merged.push(first.clone());
    }
    for block in blocks.iter().skip(1) {
        let last = merged.last_mut().expect("merged non-empty");
        if block.start <= last.end {
            last.end = last.end.max(block.end);
        } else {
            merged.push(block.clone());
        }
    }
    merged
}

/// Split a .loc index into chunks of approximately `chunk_size` bytes.
pub fn split_loc_file(
    loc_file: &str,
    chunk_size: usize,
) -> anyhow::Result<Vec<(String, u64, usize)>> {
    let loc_of = load_loc(loc_file)?;

    let mut chunks: Vec<(String, u64, usize)> = Vec::new();
    let mut cur_size = 0;
    let mut cur_start_offset = 0;
    let mut cur_first_seq = String::new();

    for (seq_id, &(offset, size)) in &loc_of {
        if cur_size + size > chunk_size && !cur_first_seq.is_empty() {
            chunks.push((cur_first_seq.clone(), cur_start_offset, cur_size));
            cur_size = 0;
            cur_start_offset = offset;
            cur_first_seq = seq_id.clone();
        }

        if cur_size == 0 {
            cur_start_offset = offset;
            cur_first_seq = seq_id.clone();
        }

        cur_size += size;
    }

    if !cur_first_seq.is_empty() {
        chunks.push((cur_first_seq, cur_start_offset, cur_size));
    }

    Ok(chunks)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{Duration, SystemTime};

    /// Regression: a `.loc` index older than its FASTA was served as-is, so
    /// an edited genome returned stale offsets/sizes (slice errors, or worse,
    /// wrong sequence on same-length edits). The index must be rebuilt.
    #[test]
    fn stale_loc_index_is_rebuilt() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("g.fa");
        let loc = dir.path().join("g.fa.loc");
        std::fs::write(&fa, format!(">chr1\n{}\n", "ACGT".repeat(250))).unwrap();

        let (mut reader, loc_of) = open_indexed(fa.to_str().unwrap(), false).unwrap();
        let seq = fetch_range_seq(&mut reader, &loc_of, &Range::from("chr1", 1, 100)).unwrap();
        assert_eq!(seq, "ACGT".repeat(25));

        // Edit the FASTA in place (same length, different content) and age
        // the index so its mtime precedes the FASTA's.
        std::fs::write(&fa, format!(">chr1\n{}\n", "TGCA".repeat(250))).unwrap();
        std::fs::File::open(&loc)
            .unwrap()
            .set_modified(SystemTime::now() - Duration::from_secs(10))
            .unwrap();

        let (mut reader, loc_of) = open_indexed(fa.to_str().unwrap(), false).unwrap();
        let seq = fetch_range_seq(&mut reader, &loc_of, &Range::from("chr1", 1, 100)).unwrap();
        assert_eq!(seq, "TGCA".repeat(25));
    }
}
