use std::io::{self, Write};

/// A single PAF (Pairwise mApping Format) record with 12 mandatory columns
/// and optional SAM-like tags.
///
/// Follows the [lh3/miniasm PAF specification](https://github.com/lh3/miniasm/blob/master/PAF.md).
#[derive(Debug, Clone, PartialEq)]
pub struct PafRecord {
    /// Column 1: Query sequence name.
    pub query_name: String,
    /// Column 2: Query sequence length.
    pub query_length: u32,
    /// Column 3: Query start (0-based, inclusive).
    pub query_start: u32,
    /// Column 4: Query end (0-based, exclusive).
    pub query_end: u32,
    /// Column 5: Strand (`+` forward, `-` reverse complement).
    pub strand: char,
    /// Column 6: Target sequence name.
    pub target_name: String,
    /// Column 7: Target sequence length.
    pub target_length: u32,
    /// Column 8: Target start (0-based, inclusive).
    pub target_start: u32,
    /// Column 9: Target end (0-based, exclusive).
    pub target_end: u32,
    /// Column 10: Number of matching bases.
    pub matches: u32,
    /// Column 11: Alignment block length.
    pub block_length: u32,
    /// Column 12: Mapping quality (255 if unavailable).
    pub mapq: u8,
    /// Optional SAM-like tags (e.g. `gi:f:0.95`, `cg:Z:10M`).
    pub tags: Vec<String>,
}

/// Write a PAF record with all tags.
///
/// Tags are written tab-separated after the 12 mandatory columns.
pub fn write_paf_record<W: Write>(writer: &mut W, rec: &PafRecord) -> io::Result<()> {
    write!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        rec.query_name,
        rec.query_length,
        rec.query_start,
        rec.query_end,
        rec.strand,
        rec.target_name,
        rec.target_length,
        rec.target_start,
        rec.target_end,
        rec.matches,
        rec.block_length,
        rec.mapq,
    )?;
    for tag in &rec.tags {
        write!(writer, "\t{}", tag)?;
    }
    writeln!(writer)
}
