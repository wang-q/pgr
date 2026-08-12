//! Quality trimming for FASTQ reads.

use std::io::Write;

use anyhow::{bail, Context};

use crate::libs::fmt::seq::{SeqReader, SeqRecord};
use crate::libs::fq::qual;

/// FASTQ quality encoding.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QualityBase {
    /// Auto-detect `PHRED33` vs `PHRED64` from a record sample.
    Auto,
    /// Sanger / Illumina 1.8+.
    Phred33,
    /// Illumina 1.3-1.7 / Solexa.
    Phred64,
}

impl QualityBase {
    /// Returns the numeric Phred offset, auto-detecting from `sample` when `Auto`.
    pub fn offset(self, sample: &[SeqRecord]) -> u8 {
        match self {
            QualityBase::Auto => qual::detect_quality_base(sample),
            QualityBase::Phred33 => qual::PHRED33,
            QualityBase::Phred64 => qual::PHRED64,
        }
    }
}

/// Quality trimming algorithm.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Method {
    /// Sickle-style adaptive sliding window.
    Sliding,
    /// cutadapt-style Mott cumulative quality.
    Mott,
}

/// Options for `fq trim-qual`.
#[derive(Clone, Copy, Debug)]
pub struct TrimOptions {
    /// Quality threshold (Phred score).
    pub qual_threshold: f64,
    /// Minimum kept length; shorter reads are discarded.
    pub length_threshold: usize,
    /// Trimming algorithm.
    pub method: Method,
    /// Disable 5' trimming.
    pub no_fiveprime: bool,
    /// Input quality encoding.
    pub quality_base: QualityBase,
    /// Trim 3' poly-G tails of at least this length (0 disables).
    pub polyg_right: usize,
}

/// Records sampled for quality-base auto-detection.
const SAMPLE_READS: usize = 1000;
/// Quality bytes scanned before detection stops (1 MB).
const SAMPLE_QUAL_BYTES: usize = 1 << 20;

/// Reads up to `SAMPLE_READS` records (or `SAMPLE_QUAL_BYTES` quality bytes)
/// from `reader` into `sample`.
fn sample_records(reader: &mut SeqReader, sample: &mut Vec<SeqRecord>) -> anyhow::Result<()> {
    let mut total_qual = 0usize;
    let mut rec = SeqRecord::new();
    while sample.len() < SAMPLE_READS && total_qual < SAMPLE_QUAL_BYTES {
        if !reader.read_record(&mut rec)? {
            break;
        }
        total_qual += rec.quality_scores().len();
        sample.push(rec.clone());
    }
    Ok(())
}

/// Sickle-style adaptive sliding-window cuts (`sliding.c` `sliding_window`).
///
/// Window size is 10% of the read length (min 1). The 5' cut is the first
/// base reaching the threshold in the first window whose average is >= it;
/// the 3' cut is the first base below the threshold in the next window whose
/// average drops below it. Returns `None` when no 5' cut is found.
fn sliding_cut(
    qual: &[u8],
    base: u8,
    threshold: f64,
    no_fiveprime: bool,
) -> Option<(usize, usize)> {
    let n = qual.len();
    if n == 0 {
        return Some((0, 0));
    }
    let window_size = usize::max(1, n / 10);
    let q = |i: usize| qual[i] as f64 - base as f64;
    let mut window_total: f64 = (0..window_size).map(&q).sum();
    let mut window_start = 0usize;
    let mut five = 0usize;
    let mut three = n;
    let mut found_five = false;
    while window_start + window_size <= n {
        let avg = window_total / window_size as f64;
        if !no_fiveprime && !found_five && avg >= threshold {
            five = (window_start..window_start + window_size)
                .find(|&j| q(j) >= threshold)
                .unwrap_or(window_start);
            found_five = true;
        }
        if avg < threshold && (found_five || no_fiveprime) {
            three = (window_start..window_start + window_size)
                .find(|&j| q(j) < threshold)
                .unwrap_or(three);
            break;
        }
        window_total -= q(window_start);
        if window_start + window_size < n {
            window_total += q(window_start + window_size);
        }
        window_start += 1;
    }
    if !found_five && !no_fiveprime {
        None
    } else {
        Some((five, three))
    }
}

/// cutadapt-style Mott cuts (`qualtrim.pyx` `quality_trim_index`).
///
/// Independent 5'/3' scans accumulate `cutoff - quality`; the cut is the point
/// of maximal cumulative score before the sum turns negative. Returns `(0, 0)`
/// for an invalid (empty) interval.
fn mott_cut(qual: &[u8], base: u8, cutoff_front: f64, cutoff_back: f64) -> (usize, usize) {
    let n = qual.len();
    let score = |q: u8| cutoff_front - (q as f64 - base as f64);
    let (mut start, mut s, mut max_s) = (0usize, 0f64, 0f64);
    for (i, &qc) in qual.iter().enumerate() {
        s += score(qc);
        if s < 0.0 {
            break;
        }
        if s > max_s {
            max_s = s;
            start = i + 1;
        }
    }
    let score = |q: u8| cutoff_back - (q as f64 - base as f64);
    let (mut stop, mut s, mut max_s) = (n, 0f64, 0f64);
    for (i, &qc) in qual.iter().enumerate().rev() {
        s += score(qc);
        if s < 0.0 {
            break;
        }
        if s > max_s {
            max_s = s;
            stop = i;
        }
    }
    if start >= stop {
        (0, 0)
    } else {
        (start, stop)
    }
}

/// Trims a 3' poly-G tail of at least `min_run` consecutive Gs from `end`
/// (BBDuk `trimpolygright` with `maxNonPoly=0`).
fn polyg_end(seq: &[u8], end: usize, min_run: usize) -> usize {
    if min_run == 0 {
        return end;
    }
    let mut run = 0;
    let mut i = end;
    while i > 0 && seq[i - 1] == b'G' {
        run += 1;
        i -= 1;
    }
    if run >= min_run {
        i
    } else {
        end
    }
}

/// Checks that every quality character decodes to a Phred score in [0, 93].
fn validate_quality(rec: &SeqRecord, base: u8) -> anyhow::Result<()> {
    for (i, &qc) in rec.quality_scores().iter().enumerate() {
        let q = qc as i32 - base as i32;
        if !(0..=93).contains(&q) {
            bail!(
                "invalid quality character {:?} at position {} for quality base {} in record {:?}",
                qc as char,
                i,
                base,
                rec.name()
            );
        }
    }
    Ok(())
}

/// Computes the kept interval `[start, end)` for a record, or `None` when the
/// read must be discarded (too short, or no valid interval).
fn trim_interval(
    rec: &SeqRecord,
    base: u8,
    opts: &TrimOptions,
) -> anyhow::Result<Option<(usize, usize)>> {
    let seq = rec.sequence();
    let qual = rec.quality_scores();
    validate_quality(rec, base)?;
    if seq.len() < opts.length_threshold {
        return Ok(None);
    }
    let (start, end) = match opts.method {
        Method::Sliding => match sliding_cut(qual, base, opts.qual_threshold, opts.no_fiveprime) {
            Some(cut) => cut,
            None => return Ok(None),
        },
        Method::Mott => {
            let front = if opts.no_fiveprime {
                0.0
            } else {
                opts.qual_threshold
            };
            mott_cut(qual, base, front, opts.qual_threshold)
        }
    };
    let end = polyg_end(seq, end, opts.polyg_right);
    if end < start || end - start < opts.length_threshold {
        return Ok(None);
    }
    Ok(Some((start, end)))
}

/// Writes a FASTQ record keeping the original name and comment.
fn write_record<W: Write>(
    w: &mut W,
    rec: &SeqRecord,
    seq: &[u8],
    qual: &[u8],
) -> std::io::Result<()> {
    let comment = rec.comment();
    if comment.is_empty() {
        writeln!(w, "@{}", rec.name())?;
    } else {
        writeln!(w, "@{} {}", rec.name(), comment)?;
    }
    w.write_all(seq)?;
    w.write_all(b"\n+\n")?;
    w.write_all(qual)?;
    w.write_all(b"\n")?;
    Ok(())
}

/// Trims one record and writes it when it survives.
fn process_record<W: Write>(
    rec: &SeqRecord,
    base: u8,
    opts: &TrimOptions,
    out: &mut W,
) -> anyhow::Result<()> {
    if !rec.is_fastq() {
        bail!("input is not FASTQ (record {:?})", rec.name());
    }
    if let Some((start, end)) = trim_interval(rec, base, opts)? {
        write_record(
            out,
            rec,
            &rec.sequence()[start..end],
            &rec.quality_scores()[start..end],
        )?;
    }
    Ok(())
}

/// Trims one read pair and routes the survivors: both kept go to `out1` (and
/// `out2` when given, otherwise interleaved into `out1`); a single survivor
/// goes to `single` when given.
#[allow(clippy::too_many_arguments)]
fn process_pair<W: Write>(
    r1: &SeqRecord,
    r2: &SeqRecord,
    base: u8,
    opts: &TrimOptions,
    out1: &mut W,
    out2: Option<&mut &mut W>,
    single: Option<&mut &mut W>,
) -> anyhow::Result<()> {
    let keep1 = trim_interval(r1, base, opts)?;
    let keep2 = trim_interval(r2, base, opts)?;
    match (keep1, keep2) {
        (Some((s1, e1)), Some((s2, e2))) => {
            write_record(
                out1,
                r1,
                &r1.sequence()[s1..e1],
                &r1.quality_scores()[s1..e1],
            )?;
            let seq2 = &r2.sequence()[s2..e2];
            let qual2 = &r2.quality_scores()[s2..e2];
            match out2 {
                Some(w) => write_record(w, r2, seq2, qual2)?,
                None => write_record(out1, r2, seq2, qual2)?,
            }
        }
        (Some((s1, e1)), None) => {
            if let Some(w) = single {
                write_record(w, r1, &r1.sequence()[s1..e1], &r1.quality_scores()[s1..e1])?;
            }
        }
        (None, Some((s2, e2))) => {
            if let Some(w) = single {
                write_record(w, r2, &r2.sequence()[s2..e2], &r2.quality_scores()[s2..e2])?;
            }
        }
        (None, None) => {}
    }
    Ok(())
}

/// Runs quality trimming on a single-end FASTQ file.
pub fn run_single<W: Write>(infile: &str, out: &mut W, opts: &TrimOptions) -> anyhow::Result<()> {
    let mut reader =
        SeqReader::new(infile).with_context(|| format!("Failed to open reader for {}", infile))?;
    let mut sample = Vec::new();
    if opts.quality_base == QualityBase::Auto {
        sample_records(&mut reader, &mut sample)?;
    }
    let base = opts.quality_base.offset(&sample);
    for rec in &sample {
        process_record(rec, base, opts, out)?;
    }
    let mut rec = SeqRecord::new();
    while reader.read_record(&mut rec)? {
        process_record(&rec, base, opts, out)?;
    }
    Ok(())
}

/// Runs quality trimming on paired-end FASTQ files.
///
/// Quality encoding is detected from `infile1` only. When the files have
/// different record counts, a warning is printed and only the common prefix
/// is processed.
pub fn run_paired<W: Write>(
    infile1: &str,
    infile2: &str,
    out1: &mut W,
    mut out2: Option<&mut W>,
    mut single: Option<&mut W>,
    opts: &TrimOptions,
) -> anyhow::Result<()> {
    let mut reader1 = SeqReader::new(infile1)
        .with_context(|| format!("Failed to open reader for {}", infile1))?;
    let mut reader2 = SeqReader::new(infile2)
        .with_context(|| format!("Failed to open reader for {}", infile2))?;
    let mut sample = Vec::new();
    if opts.quality_base == QualityBase::Auto {
        sample_records(&mut reader1, &mut sample)?;
    }
    let base = opts.quality_base.offset(&sample);

    let mut warned = false;
    let mut rec1 = SeqRecord::new();
    let mut rec2 = SeqRecord::new();
    for r1 in &sample {
        if !reader2.read_record(&mut rec2)? {
            warn_pair_mismatch(&mut warned);
            return Ok(());
        }
        process_pair(r1, &rec2, base, opts, out1, out2.as_mut(), single.as_mut())?;
    }
    loop {
        let has1 = reader1.read_record(&mut rec1)?;
        let has2 = reader2.read_record(&mut rec2)?;
        if !has1 && !has2 {
            break;
        }
        if !has1 || !has2 {
            warn_pair_mismatch(&mut warned);
            break;
        }
        process_pair(
            &rec1,
            &rec2,
            base,
            opts,
            out1,
            out2.as_mut(),
            single.as_mut(),
        )?;
    }
    Ok(())
}

/// Prints the paired-count mismatch warning once.
fn warn_pair_mismatch(warned: &mut bool) {
    if !*warned {
        eprintln!("warning: paired input files have different numbers of reads; processing common prefix only");
        *warned = true;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn records_from(fastq: &str) -> Vec<SeqRecord> {
        let mut reader = SeqReader::from_reader(Box::new(Cursor::new(fastq.as_bytes())));
        let mut rec = SeqRecord::new();
        let mut out = Vec::new();
        while reader.read_record(&mut rec).unwrap() {
            out.push(rec.clone());
        }
        out
    }

    #[test]
    fn sliding_keeps_high_quality_read() {
        // All Q30 ('?'), length 40; window 4, average never below 20.
        let rec = &records_from("@r\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\n????????????????????????????????????????\n")[0];
        assert_eq!(
            sliding_cut(rec.quality_scores(), 33, 20.0, false),
            Some((0, 40))
        );
    }

    #[test]
    fn sliding_trims_low_quality_three_prime() {
        // 30 good bases then 10 bad ones; window 4; 3' cut inside the bad run.
        let seq = "A".repeat(40);
        let qual: String = "?".repeat(30) + &"!".repeat(10);
        let fastq = format!("@r\n{seq}\n+\n{qual}\n");
        let rec = &records_from(&fastq)[0];
        let (five, three) = sliding_cut(rec.quality_scores(), 33, 20.0, false).unwrap();
        assert_eq!(five, 0);
        assert_eq!(three, 30);
    }

    #[test]
    fn sliding_discards_when_no_five_prime_found() {
        // All Q10 ('+'); window average never reaches 20.
        let rec = &records_from("@r\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\n++++++++++++++++++++++++++++++++++++++++\n")[0];
        assert_eq!(sliding_cut(rec.quality_scores(), 33, 20.0, false), None);
    }

    #[test]
    fn mott_trims_low_quality_three_prime() {
        // 30 Q30 bases then 10 Q10 bases; stop ends before the bad tail.
        let seq = "A".repeat(40);
        let qual: String = "?".repeat(30) + &"+".repeat(10);
        let rec = &records_from(&format!("@r\n{seq}\n+\n{qual}\n"))[0];
        let (start, stop) = mott_cut(rec.quality_scores(), 33, 20.0, 20.0);
        assert_eq!(start, 0);
        assert!(stop <= 30);
    }

    #[test]
    fn mott_keeps_high_quality_read() {
        let rec = &records_from("@r\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\n????????????????????????????????????????\n")[0];
        assert_eq!(mott_cut(rec.quality_scores(), 33, 20.0, 20.0), (0, 40));
    }

    #[test]
    fn mott_discards_all_low_quality_read() {
        let rec = &records_from("@r\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\n++++++++++++++++++++++++++++++++++++++++\n")[0];
        assert_eq!(mott_cut(rec.quality_scores(), 33, 20.0, 20.0), (0, 0));
    }

    #[test]
    fn polyg_trims_only_long_enough_tails() {
        assert_eq!(polyg_end(b"ACGGGG", 6, 3), 2);
        assert_eq!(polyg_end(b"ACGGGG", 6, 7), 6);
        assert_eq!(polyg_end(b"ACGGGG", 6, 0), 6);
    }

    #[test]
    fn interval_drops_short_reads() {
        let opts = TrimOptions {
            qual_threshold: 20.0,
            length_threshold: 20,
            method: Method::Sliding,
            no_fiveprime: false,
            quality_base: QualityBase::Phred33,
            polyg_right: 0,
        };
        let rec = &records_from("@r\nACGT\n+\n!!!!\n")[0];
        assert!(trim_interval(rec, 33, &opts).unwrap().is_none());
    }

    #[test]
    fn invalid_quality_errors() {
        // ASCII 32 (space) decodes below 0 for base 33.
        let rec = &records_from("@r\nACGT\n+\n! !!\n")[0];
        assert!(trim_interval(
            rec,
            33,
            &TrimOptions {
                qual_threshold: 20.0,
                length_threshold: 1,
                method: Method::Sliding,
                no_fiveprime: false,
                quality_base: QualityBase::Phred33,
                polyg_right: 0,
            }
        )
        .is_err());
    }
}
