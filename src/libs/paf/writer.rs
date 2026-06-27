use super::record::PafRecord;
use std::io::{self, Write};

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
