//! FASTQ record writers shared by `fq` commands.

use std::io::{self, Read, Write};

use anyhow::{bail, Context};

/// Detect whether `path` is a FASTQ file (as opposed to FASTA) by inspecting
/// the first byte of the content (after gzip decompression if needed).
pub fn is_fq<P: AsRef<std::path::Path>>(path: P) -> anyhow::Result<bool> {
    let path = path.as_ref();

    let mut buffer = [0; 2];
    {
        let mut file = std::fs::File::open(path)
            .with_context(|| format!("could not open {}", path.display()))?;
        file.read_exact(&mut buffer)
            .context("could not read header")?;
    }

    let first_char = if buffer[0] == 0x1f && buffer[1] == 0x8b {
        // Gzip-compressed: decompress the first bytes
        let mut decoder = flate2::read::GzDecoder::new(
            std::fs::File::open(path)
                .with_context(|| format!("could not open {}", path.display()))?,
        );
        let mut buf = [0; 2];
        decoder
            .read_exact(&mut buf)
            .context("could not read decompressed header")?;
        buf[0] as char
    } else {
        buffer[0] as char
    };

    match first_char {
        '>' => Ok(false),
        '@' => Ok(true),
        c => bail!("unknown file format, leading byte: {:?}", c),
    }
}

/// Write a FASTQ record (4-line form) to `writer`.
pub fn write_fq<W: Write>(writer: &mut W, name: &str, seq: &[u8], qual: &[u8]) -> io::Result<()> {
    writer.write_fmt(format_args!("@{}\n", name))?;
    writer.write_all(seq)?;
    writer.write_all(b"\n")?;
    writer.write_all(b"+\n")?;
    writer.write_all(qual)?;
    writer.write_all(b"\n")?;
    Ok(())
}

/// Write a FASTA record (2-line form) to `writer`.
pub fn write_fa<W: Write>(writer: &mut W, name: &str, seq: &[u8]) -> io::Result<()> {
    writer.write_fmt(format_args!(">{}\n", name))?;
    writer.write_all(seq)?;
    writer.write_all(b"\n")?;
    Ok(())
}

/// Write a paired-end record pair (R1 then R2) in either FASTQ or FASTA form.
///
/// `*_qual` is `Some` when the source record carries quality (FASTQ input);
/// for FASTQ output with `None`, quality is filled with `'!'` of the same
/// length as the corresponding sequence. For FASTA output, quality is
/// ignored.
#[allow(clippy::too_many_arguments)]
pub fn write_pair<W: Write>(
    writer: &mut W,
    prefix: &str,
    idx: usize,
    r1_seq: &[u8],
    r1_qual: Option<&[u8]>,
    r2_seq: &[u8],
    r2_qual: Option<&[u8]>,
    is_out_fq: bool,
) -> anyhow::Result<()> {
    let r1_name = format!("{}{}/1", prefix, idx);
    let r2_name = format!("{}{}/2", prefix, idx);
    if is_out_fq {
        let r1_default: Vec<u8> = vec![b'!'; r1_seq.len()];
        let r1_q = r1_qual.unwrap_or(&r1_default);
        let r2_default: Vec<u8> = vec![b'!'; r2_seq.len()];
        let r2_q = r2_qual.unwrap_or(&r2_default);
        write_fq(writer, &r1_name, r1_seq, r1_q)?;
        write_fq(writer, &r2_name, r2_seq, r2_q)?;
    } else {
        write_fa(writer, &r1_name, r1_seq)?;
        write_fa(writer, &r2_name, r2_seq)?;
    }
    Ok(())
}

/// Interleave paired-end reads. Single-file mode generates dummy R2.
/// Returns the final read index (start + count).
pub fn interleave<W: Write>(
    writer: &mut W,
    infiles: &[String],
    prefix: &str,
    start: usize,
    is_out_fq: bool,
) -> anyhow::Result<usize> {
    let is_in_fq = is_fq(&infiles[0])?;
    let mut idx = start;

    if infiles.len() == 1 {
        let infile = &infiles[0];
        if is_in_fq {
            let reader = crate::libs::io::reader(infile)?;
            let mut seq_in = noodles_fastq::io::Reader::new(reader);
            for result in seq_in.records() {
                let record = result?;
                // Dummy R2 is always a single "N"; quality is ignored for FA
                // output and "!" for FQ output (via write_pair's default).
                let r2_seq: &[u8] = b"N";
                write_pair(
                    writer,
                    prefix,
                    idx,
                    record.sequence(),
                    Some(record.quality_scores()),
                    r2_seq,
                    Some(b"!"),
                    is_out_fq,
                )?;
                idx += 1;
            }
        } else {
            let mut seq_in = crate::libs::fmt::fa::reader(infile)?;
            for result in seq_in.records() {
                let record = result?;
                write_pair(
                    writer,
                    prefix,
                    idx,
                    &record.sequence()[..],
                    None,
                    b"N",
                    None,
                    is_out_fq,
                )?;
                idx += 1;
            }
        }
    } else {
        if is_in_fq {
            let reader1 = crate::libs::io::reader(&infiles[0])?;
            let mut seq1_in = noodles_fastq::io::Reader::new(reader1);
            let reader2 = crate::libs::io::reader(&infiles[1])?;
            let mut seq2_in = noodles_fastq::io::Reader::new(reader2);
            idx = interleave_pair_iter(
                writer,
                &mut seq1_in.records(),
                &mut seq2_in.records(),
                prefix,
                idx,
                is_out_fq,
                |record1, record2| {
                    (
                        record1.sequence(),
                        Some(record1.quality_scores()),
                        record2.sequence(),
                        Some(record2.quality_scores()),
                    )
                },
            )?;
        } else {
            let mut seq1_in = crate::libs::fmt::fa::reader(&infiles[0])?;
            let mut seq2_in = crate::libs::fmt::fa::reader(&infiles[1])?;
            idx = interleave_pair_iter(
                writer,
                &mut seq1_in.records(),
                &mut seq2_in.records(),
                prefix,
                idx,
                is_out_fq,
                |record1, record2| (&record1.sequence()[..], None, &record2.sequence()[..], None),
            )?;
        }
    }

    Ok(idx)
}

/// Interleave two record streams pair-by-pair, erroring if the two files hold
/// a different number of reads instead of silently truncating to the shorter.
fn interleave_pair_iter<W, I1, I2, R1, R2, E, F>(
    writer: &mut W,
    iter1: &mut I1,
    iter2: &mut I2,
    prefix: &str,
    mut idx: usize,
    is_out_fq: bool,
    map: F,
) -> anyhow::Result<usize>
where
    W: Write,
    I1: Iterator<Item = Result<R1, E>>,
    I2: Iterator<Item = Result<R2, E>>,
    E: std::error::Error + Send + Sync + 'static,
    F: for<'a> Fn(&'a R1, &'a R2) -> (&'a [u8], Option<&'a [u8]>, &'a [u8], Option<&'a [u8]>),
{
    loop {
        let r1 = iter1.next();
        let r2 = iter2.next();
        match (r1, r2) {
            (Some(Ok(record1)), Some(Ok(record2))) => {
                let (seq1, qual1, seq2, qual2) = map(&record1, &record2);
                write_pair(writer, prefix, idx, seq1, qual1, seq2, qual2, is_out_fq)?;
                idx += 1;
            }
            (None, None) => break,
            (Some(Err(e)), _) | (_, Some(Err(e))) => return Err(e.into()),
            _ => {
                anyhow::bail!("paired files have different numbers of reads");
            }
        }
    }
    Ok(idx)
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn test_is_fq_plain_text() {
        let dir = tempdir().unwrap();

        // Create a plain text FASTQ file
        let fq_file_path = dir.path().join("test.fq");
        {
            let mut file = std::fs::File::create(&fq_file_path).unwrap();
            writeln!(file, "@SEQ_ID").unwrap(); // FASTQ format
        }
        assert!(is_fq(&fq_file_path).unwrap());

        // Create a plain text FASTA file
        let fasta_file_path = dir.path().join("test.fasta");
        {
            let mut file = std::fs::File::create(&fasta_file_path).unwrap();
            writeln!(file, ">SEQ_ID").unwrap(); // FASTA format
        }
        assert!(!is_fq(&fasta_file_path).unwrap());
    }

    #[test]
    fn test_is_fq_gzip() {
        let dir = tempdir().unwrap();

        // Create a Gzip FASTQ file
        let fq_file_path = dir.path().join("test.fq.gz");
        {
            let file = std::fs::File::create(&fq_file_path).unwrap();
            let mut encoder = GzEncoder::new(file, flate2::Compression::default());
            writeln!(encoder, "@SEQ_ID").unwrap(); // FASTQ format
            encoder.finish().unwrap();
        }
        assert!(is_fq(&fq_file_path).unwrap());

        // Create a Gzip FASTA file
        let fasta_file_path = dir.path().join("test.fasta.gz");
        {
            let file = std::fs::File::create(&fasta_file_path).unwrap();
            let mut encoder = GzEncoder::new(file, flate2::Compression::default());
            writeln!(encoder, ">SEQ_ID").unwrap(); // FASTA format
            encoder.finish().unwrap();
        }
        assert!(!is_fq(&fasta_file_path).unwrap());
    }

    #[test]
    fn test_interleave_two_files_returns_final_index() {
        // Regression: the two-file path used to discard the updated index from
        // interleave_pair_iter, so `interleave` returned only `start` instead
        // of the documented final read index (start + count).
        let dir = tempdir().unwrap();
        let r1 = dir.path().join("R1.fq");
        let r2 = dir.path().join("R2.fq");
        std::fs::write(&r1, "@R1\nACGT\n+\n!!!!\n@R3\nACGT\n+\n!!!!\n").unwrap();
        std::fs::write(&r2, "@R2\nACGT\n+\n!!!!\n@R4\nACGT\n+\n!!!!\n").unwrap();

        let mut out = Vec::new();
        let final_idx = interleave(
            &mut out,
            &[
                r1.to_str().unwrap().to_string(),
                r2.to_str().unwrap().to_string(),
            ],
            "read",
            5,
            false,
        )
        .unwrap();
        assert_eq!(final_idx, 7); // start(5) + 2 pairs
    }
}
