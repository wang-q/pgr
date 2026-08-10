//! Split an interleaved FASTQ stream into R1/R2/singles files.

use crate::libs::fmt::fq::write_fq;
use crate::libs::fmt::seq::{SeqReader, SeqRecord};
use anyhow::Result;
use std::io::Write;

/// Splits interleaved FASTQ `infile` into paired records (`out1`/`out2`) and
/// the trailing unpaired record (written to `single` when provided).
///
/// Headers, order, and formatting are preserved; a record without its mate is
/// sent to `single` or silently dropped when no singles writer is given.
pub fn split<W: Write>(
    infile: &str,
    out1: &mut W,
    out2: &mut W,
    mut single: Option<&mut W>,
) -> Result<()> {
    let mut reader = SeqReader::new(infile)?;
    let mut rec1 = SeqRecord::new();
    let mut rec2 = SeqRecord::new();
    loop {
        if !reader.read_record(&mut rec1)? {
            break;
        }
        if !reader.read_record(&mut rec2)? {
            if let Some(w) = single.take() {
                write_record(w, &rec1)?;
            } else {
                eprintln!("warning: unpaired read discarded (no --outfile-single)");
            }
            break;
        }
        write_record(out1, &rec1)?;
        write_record(out2, &rec2)?;
    }
    Ok(())
}

/// Writes a FASTQ record, preserving the `name comment` header layout.
fn write_record<W: Write>(w: &mut W, rec: &SeqRecord) -> anyhow::Result<()> {
    let comment = rec.comment();
    let header = if comment.is_empty() {
        rec.name().to_string()
    } else {
        format!("{} {}", rec.name(), comment)
    };
    write_fq(w, &header, rec.sequence(), rec.quality_scores())?;
    Ok(())
}
