//! SAM alignment parsing and conversion (via `noodles-sam`).
//!
//! Reads SAM through the noodles-sam streaming reader and converts mapped
//! alignments to `.rg` range lines for `pgr sam to-rg`.

use anyhow::{bail, Context, Result};
use noodles_sam::alignment::RecordBuf;
use std::io::{BufRead, Write};

/// Converts SAM alignments to `.rg` range lines (`chr:start-end`, 1-based
/// inclusive). Unmapped records are skipped; malformed records are skipped
/// unless `strict`.
pub fn to_ranges<R: BufRead, W: Write>(reader: R, writer: &mut W, strict: bool) -> Result<()> {
    let mut reader = noodles_sam::io::Reader::new(reader);
    reader.read_header().context("failed to read SAM header")?;
    let mut record = noodles_sam::Record::default();
    loop {
        match reader.read_record(&mut record) {
            Ok(0) => break, // EOF
            Ok(_) => {
                let flags = match record.flags() {
                    Ok(flags) => flags,
                    Err(e) if strict => bail!("malformed SAM FLAG: {e}"),
                    Err(_) => continue,
                };
                if flags.is_unmapped() {
                    continue;
                }
                let Some(rname) = record.reference_sequence_name() else {
                    continue;
                };
                let Some(start) = record.alignment_start() else {
                    continue;
                };
                let start = match start {
                    Ok(start) => start,
                    Err(e) if strict => bail!("malformed SAM POS: {e}"),
                    Err(_) => continue,
                };
                let mut span = 0usize;
                let mut bad = false;
                for op in record.cigar().iter() {
                    match op {
                        Ok(op) if op.kind().consumes_reference() => span += op.len(),
                        Ok(_) => {}
                        Err(_) => {
                            bad = true;
                            break;
                        }
                    }
                }
                if bad {
                    if strict {
                        bail!("malformed SAM CIGAR");
                    }
                    continue;
                }
                if span == 0 {
                    continue;
                }
                let start = start.get();
                writeln!(writer, "{}:{}-{}", rname, start, start + span - 1)?;
            }
            Err(e) if strict => bail!("malformed SAM record: {e}"),
            Err(_) => continue, // the malformed line was consumed; keep reading
        }
    }
    Ok(())
}

/// Computes the insert-size histogram of a paired SAM and writes it in the
/// BBTools `reformat.sh ihist` text format (`#Mean/#Median/#Mode/#STDev/
/// #PercentOfPairs` + `#InsertSize\tCount` rows).
///
/// Pairs are grouped by read name (first whitespace token, trailing
/// `/1`/`/2` stripped); only proper FR pairs — both ends mapped, same
/// reference, opposite strands, pointing inward — contribute an insert
/// size. `#PercentOfPairs` is the fraction of pairs contributing to the
/// histogram.
pub fn ihist<R: BufRead, W: Write>(reader: R, writer: &mut W) -> Result<()> {
    use std::collections::HashMap;

    let mut reader = noodles_sam::io::Reader::new(reader);
    let header = reader.read_header().context("failed to read SAM header")?;

    // Pending first/last segment per normalized pair name.
    let mut pending: HashMap<Vec<u8>, [Option<RecordBuf>; 2]> = HashMap::new();
    let mut total_pairs = 0u64;
    let mut insert_sizes: Vec<u32> = Vec::new();

    for result in reader.record_bufs(&header) {
        let record = result.context("failed to read SAM record")?;
        let flags = record.flags();
        if !flags.is_segmented() {
            continue;
        }
        let Some(name) = record.name() else {
            continue;
        };
        let key = pair_key(name);
        let slot = if flags.is_first_segment() { 0 } else { 1 };
        let entry = pending.entry(key.clone()).or_default();
        if entry[slot].is_none() {
            entry[slot] = Some(record);
        }
        if flags.is_first_segment() {
            total_pairs += 1;
        }
        if entry[0].is_some() && entry[1].is_some() {
            let r1 = entry[0].take().unwrap();
            let r2 = entry[1].take().unwrap();
            pending.remove(&key);
            if let Some(size) = proper_insert_size(&r1, &r2) {
                insert_sizes.push(size);
            }
        }
    }

    insert_sizes.sort_unstable();
    let n = insert_sizes.len() as u64;
    let (mean, median, mode, stdev) = if insert_sizes.is_empty() {
        (0.0, 0u64, 0u64, 0.0)
    } else {
        let sum: u64 = insert_sizes.iter().map(|&x| x as u64).sum();
        let mean = sum as f64 / n as f64;
        let median = insert_sizes[((n - 1) / 2) as usize] as u64;
        // Mode: most frequent insert size; ties -> the smallest.
        let mut mode = insert_sizes[0];
        let mut best_count = 0usize;
        let mut i = 0usize;
        while i < insert_sizes.len() {
            let mut j = i;
            while j < insert_sizes.len() && insert_sizes[j] == insert_sizes[i] {
                j += 1;
            }
            if j - i > best_count {
                best_count = j - i;
                mode = insert_sizes[i];
            }
            i = j;
        }
        let variance = insert_sizes
            .iter()
            .map(|&x| {
                let d = x as f64 - mean;
                d * d
            })
            .sum::<f64>()
            / n as f64;
        (mean, median, mode as u64, variance.sqrt())
    };

    writeln!(writer, "#Mean\t{mean:.3}")?;
    writeln!(writer, "#Median\t{median}")?;
    writeln!(writer, "#Mode\t{mode}")?;
    writeln!(writer, "#STDev\t{stdev:.3}")?;
    writeln!(
        writer,
        "#PercentOfPairs\t{:.3}",
        n as f64 / total_pairs.max(1) as f64
    )?;
    writeln!(writer, "#InsertSize\tCount")?;
    let mut i = 0usize;
    while i < insert_sizes.len() {
        let mut j = i;
        while j < insert_sizes.len() && insert_sizes[j] == insert_sizes[i] {
            j += 1;
        }
        writeln!(writer, "{}\t{}", insert_sizes[i], j - i)?;
        i = j;
    }
    Ok(())
}

/// Normalized pair key: first whitespace token, trailing `/1`/`/2` stripped.
fn pair_key(name: &[u8]) -> Vec<u8> {
    let end = name
        .iter()
        .position(|&b| b == b' ' || b == b'\t')
        .unwrap_or(name.len());
    let mut key = name[..end].to_vec();
    if key.ends_with(b"/1") || key.ends_with(b"/2") {
        key.truncate(key.len() - 2);
    }
    key
}

/// Insert size of a proper FR pair (`None` otherwise): both ends mapped,
/// same reference, opposite strands, reads pointing inward.
fn proper_insert_size(r1: &RecordBuf, r2: &RecordBuf) -> Option<u32> {
    let c1 = r1.reference_sequence_id()?;
    let c2 = r2.reference_sequence_id()?;
    if c1 != c2 {
        return None;
    }
    let p1 = r1.alignment_start()?.get() as u32;
    let p2 = r2.alignment_start()?.get() as u32;
    let rc1 = r1.flags().is_reverse_complemented();
    let rc2 = r2.flags().is_reverse_complemented();
    if rc1 == rc2 {
        return None;
    }
    let l1 = r1.sequence().len() as u32;
    let l2 = r2.sequence().len() as u32;
    let (left, right, right_len) = if !rc1 && rc2 && p1 < p2 {
        (p1, p2, l2)
    } else if rc1 && !rc2 && p2 < p1 {
        (p2, p1, l1)
    } else {
        return None; // outward orientation
    };
    Some(right + right_len - left)
}
