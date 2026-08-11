//! SAM alignment parsing and conversion.
//!
//! Minimal SAM reader for the subset needed by `pgr sam to-rg`: header
//! lines are skipped and mapped alignments become `.rg` range lines.

use anyhow::{bail, Result};
use std::io::{BufRead, Write};

/// Reference-consuming CIGAR operations (each advances the reference).
const REF_CONSUMING: [u8; 5] = *b"MDN=X";

/// Converts SAM alignments to `.rg` range lines (`chr:start-end`, 1-based
/// inclusive). Header and unmapped records are skipped; malformed lines are
/// skipped unless `strict`.
pub fn to_ranges<R: BufRead, W: Write>(reader: R, writer: &mut W, strict: bool) -> Result<()> {
    for line in reader.lines() {
        let line = line?;
        let line = line.trim_end();
        if line.is_empty() || line.starts_with('@') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let flag = match fields.get(1).and_then(|f| f.parse::<u16>().ok()) {
            Some(flag) => flag,
            None if strict => bail!("malformed SAM FLAG field: {line}"),
            None => continue,
        };
        if flag & 0x4 != 0 {
            continue; // unmapped
        }
        if fields.len() < 6 || fields[2] == "*" || fields[3] == "0" {
            continue;
        }
        let pos = match fields[3].parse::<u32>() {
            Ok(pos) => pos,
            Err(_) if strict => bail!("malformed SAM POS field: {line}"),
            Err(_) => continue,
        };
        let Some(span) = cigar_span(fields[5]) else {
            if strict {
                bail!("malformed SAM CIGAR: {line}");
            }
            continue;
        };
        if span == 0 {
            continue;
        }
        writeln!(writer, "{}:{}-{}", fields[2], pos, pos + span - 1)?;
    }
    Ok(())
}

/// Total reference span of a CIGAR (sum of M/D/N/=/X lengths); `None` when
/// the CIGAR is `*` or malformed.
fn cigar_span(cigar: &str) -> Option<u32> {
    if cigar == "*" || cigar.is_empty() {
        return None;
    }
    let mut span = 0u32;
    let mut len = 0u32;
    for b in cigar.bytes() {
        if b.is_ascii_digit() {
            len = len.checked_mul(10)?.checked_add((b - b'0') as u32)?;
        } else if REF_CONSUMING.contains(&b) {
            if len == 0 {
                return None;
            }
            span = span.checked_add(len)?;
            len = 0;
        } else if matches!(b, b'I' | b'S' | b'H' | b'P') {
            if len == 0 {
                return None;
            }
            len = 0;
        } else {
            return None;
        }
    }
    (len == 0).then_some(span)
}
