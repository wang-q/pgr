//! FASTA helpers wrapping noodles_fasta with pgr I/O.

use crate::libs::io;
use noodles_fasta as fasta;

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
