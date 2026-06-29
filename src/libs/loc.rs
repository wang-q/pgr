use indexmap::IndexMap;
use noodles_bgzf as bgzf;
use noodles_core;
use noodles_fasta as fasta;
use std::io::{Read, Seek, SeekFrom};

/// Random-access reader for indexed FASTA files (plain or BGZF-compressed).
pub enum Input {
    File(std::fs::File),
    Bgzf(bgzf::io::IndexedReader<std::fs::File>),
}

pub fn create_loc(infile: &str, locfile: &str, is_bgzf: bool) -> anyhow::Result<()> {
    let mut reader: Box<dyn std::io::BufRead> = if is_bgzf {
        // http://www.htslib.org/doc/bgzip.html
        // Bgzip will attempt to ensure BGZF blocks end on a newline when the input is a text file.
        // The exception to this is where a single line is larger than a BGZF block (64Kb).
        Box::new(bgzf::io::indexed_reader::Builder::default().build_from_path(infile)?)
    } else {
        crate::libs::io::reader(infile)?
    };

    let mut writer: Box<dyn std::io::Write> =
        Box::new(std::io::BufWriter::new(std::fs::File::create(locfile)?));

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

/// Open a FASTA file with .loc index for random access.
/// Creates the .loc index if it doesn't exist (or if `force_update` is true).
/// Returns the Input reader and the loaded .loc index.
#[allow(clippy::type_complexity)]
pub fn open_indexed(
    infile: &str,
    force_update: bool,
) -> anyhow::Result<(Input, IndexMap<String, (u64, usize)>)> {
    let is_bgzf = crate::is_bgzf(infile);
    let loc_file = format!("{}.loc", infile);
    if !std::path::Path::new(&loc_file).is_file() || force_update {
        create_loc(infile, &loc_file, is_bgzf)?;
    }
    let loc_of = load_loc(&loc_file)?;
    let reader = if is_bgzf {
        Input::Bgzf(bgzf::io::indexed_reader::Builder::default().build_from_path(infile)?)
    } else {
        Input::File(std::fs::File::open(std::path::Path::new(infile))?)
    };
    Ok((reader, loc_of))
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
) -> anyhow::Result<fasta::Record> {
    let (offset, size) = loc_of
        .get(name)
        .ok_or_else(|| anyhow::anyhow!("{} not found in the .loc index file", name))?;

    let data_buf = read_offset(reader, *offset, *size)?;
    let mut fa_in = fasta::io::Reader::new(&data_buf[..]);

    fa_in.read_definition(&mut String::new())?;
    let mut buf = Vec::new();
    fa_in.read_sequence(&mut buf)?;

    let definition = fasta::record::Definition::new(name, None);
    let sequence = fasta::record::Sequence::from(buf);
    let record = fasta::Record::new(definition, sequence);

    Ok(record)
}

pub fn records_offset(
    reader: &mut Input,
    offset: u64,
    size: usize,
) -> anyhow::Result<Vec<fasta::Record>> {
    let mut records = Vec::new();

    let data_buf = read_offset(reader, offset, size)?;
    let mut fa_in = fasta::io::Reader::new(&data_buf[..]);

    for result in fa_in.records() {
        // obtain record or fail with error
        let record = result?;
        records.push(record);
    }

    Ok(records)
}

pub fn fetch_range_seq(
    reader: &mut Input,
    loc_of: &IndexMap<String, (u64, usize)>,
    rg: &intspan::Range,
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
        let seq = record
            .sequence()
            .as_ref()
            .iter()
            .map(|&b| b as char)
            .collect();
        return Ok(seq);
    }

    // slice here is 1-based
    let start = noodles_core::Position::new(*rg.start() as usize)
        .ok_or_else(|| anyhow::anyhow!("invalid start position: {}", *rg.start()))?;
    let end = noodles_core::Position::new(*rg.end() as usize)
        .ok_or_else(|| anyhow::anyhow!("invalid end position: {}", *rg.end()))?;

    let mut slice = record
        .sequence()
        .slice(start..=end)
        .ok_or_else(|| anyhow::anyhow!("slice error for [{}]", rg))?;
    if rg.strand() == "-" {
        slice = slice.complement().rev().collect::<Result<_, _>>()?;
    }

    let seq = slice.as_ref().iter().map(|&b| b as char).collect();
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
