//! FASTQ record writers shared by `fq` commands.

use std::io::{self, Write};

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
