use super::record::PafRecord;
use std::io::{self, Write};

/// Write a PAF record with standard tags (`gi`, `bi`, `cg`, optional `ms`).
pub fn write_paf_record<W: Write>(
    writer: &mut W,
    rec: &PafRecord,
    gi: f64,
    bi: f64,
    cigar: &str,
    score: Option<u64>,
) -> io::Result<()> {
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
    write!(writer, "\tgi:f:{:.6}\tbi:f:{:.6}\tcg:Z:{}", gi, bi, cigar)?;
    if let Some(s) = score {
        write!(writer, "\tms:i:{}", s)?;
    }
    writeln!(writer)
}
