//! Split an interleaved FASTQ stream into R1/R2/singles files.

use crate::libs::fmt::fq::write_fq;
use crate::libs::fmt::seq::{SeqReader, SeqRecord};
use anyhow::{Context, Result};
use std::collections::HashMap;
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

/// Split by read-name pairing (`repair.sh` `rp` mode): buffers unpaired
/// reads keyed by name prefix, matches mates as they arrive, and emits the
/// leftovers as singles.
///
/// Mirrors `jgi.SplitPairsAndSingles.repair()`: the prefix is the first
/// whitespace token (or the part before `/` for a bare `name/1` name);
/// the pair number comes from a `/1` `/2` or `1:` `2:` suffix, falling back
/// to any `/1`/`/2` in the name (which also redefines the prefix to the
/// whole part before the first `/`). A read whose mate never arrives is
/// written to `single` in insertion order. Buffers unpaired reads in
/// memory, like the reference implementation.
pub fn split_repair<W: Write>(
    infile: &str,
    out1: &mut W,
    out2: &mut W,
    mut single: Option<&mut W>,
) -> Result<()> {
    let mut reader = SeqReader::new(infile).with_context(|| format!("failed to open {infile}"))?;
    let mut rec = SeqRecord::new();
    let mut pending: HashMap<String, (SeqRecord, Option<bool>)> = HashMap::new();
    let mut order: Vec<String> = Vec::new();
    while reader.read_record(&mut rec)? {
        let (prefix, is_r2) = pair_key(&rec);
        if let Some((old, old_is_r2)) = pending.remove(&prefix) {
            // Mate found: the cached read stays R1 unless it was R2.
            let (r1, r2) = if old_is_r2 == Some(true) {
                (rec.clone(), old)
            } else {
                (old, rec.clone())
            };
            write_record(out1, &r1)?;
            write_record(out2, &r2)?;
        } else {
            pending.insert(prefix.clone(), (rec.clone(), is_r2));
            order.push(prefix);
        }
    }
    if let Some(w) = single.as_mut() {
        for key in &order {
            if let Some((r, _)) = pending.remove(key) {
                write_record(w, &r)?;
            }
        }
    } else if !pending.is_empty() {
        eprintln!(
            "warning: {} unpaired reads discarded (no --outfile-single)",
            pending.len()
        );
    }
    Ok(())
}

/// Name prefix and pair number of a read (`SplitPairsAndSingles.repair`).
fn pair_key(rec: &SeqRecord) -> (String, Option<bool>) {
    let comment = rec.comment();
    let id = if comment.is_empty() {
        rec.name().to_string()
    } else {
        format!("{} {}", rec.name(), comment)
    };
    let id = id.strip_prefix('@').unwrap_or(&id);
    let mut tokens: Vec<&str> = id.split_whitespace().collect();
    if tokens.len() == 1 && id.contains('/') {
        let slash = id.find('/').unwrap();
        tokens = vec![&id[..slash], &id[slash..]];
    }
    let mut prefix = tokens[0].to_string();
    let suffix = tokens.last().copied();
    let is_r2 = match suffix {
        Some(s) if s.starts_with("/1") || s.starts_with("1:") => Some(false),
        Some(s) if s.starts_with("/2") || s.starts_with("2:") => Some(true),
        _ if id.contains("/1") || id.contains("/2") => {
            let slash = id.find('/').unwrap();
            prefix = id[..slash].to_string();
            let after = &id[slash + 1..];
            if after.starts_with('1') {
                Some(false)
            } else if after.starts_with('2') {
                Some(true)
            } else {
                None
            }
        }
        _ => None,
    };
    (prefix, is_r2)
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
