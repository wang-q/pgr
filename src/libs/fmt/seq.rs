//! Self-contained FAFQ (FASTA/FASTQ) sequential reader.
//!
//! kseq-style design (seqtk `kseq.h`, notes/design/seq-reader.md): a single
//! reader auto-detects `>`/`@` headers, record buffers are reused across
//! reads, and line scanning uses `memchr` batches instead of per-record
//! reallocation. Supports gzip/stdin through `crate::libs::io::reader`.

use anyhow::Context;
use bstr::{BStr, BString};
use memchr::{memchr, memchr2};
use std::io::BufRead;

/// A FAFQ record with reused buffers (capacity kept across reads).
#[derive(Clone, Default)]
pub struct SeqRecord {
    name: BString,
    comment: BString,
    seq: Vec<u8>,
    qual: Vec<u8>,
    is_fastq: bool,
}

impl SeqRecord {
    /// Creates an empty record.
    pub fn new() -> Self {
        Self::default()
    }

    /// Clears all fields, keeping buffer capacities.
    fn clear(&mut self) {
        self.name.clear();
        self.comment.clear();
        self.seq.clear();
        self.qual.clear();
        self.is_fastq = false;
    }

    /// Record name (without the `>`/`@` prefix), as a byte string.
    pub fn name(&self) -> &BStr {
        BStr::new(&self.name)
    }

    /// Optional description (text after the name, before end of line).
    pub fn comment(&self) -> &BStr {
        BStr::new(&self.comment)
    }

    /// Optional description, `None` when absent (noodles-compatible).
    pub fn description(&self) -> Option<&[u8]> {
        (!self.comment.is_empty()).then_some(self.comment.as_slice())
    }

    /// Sequence bases.
    pub fn sequence(&self) -> &[u8] {
        &self.seq
    }

    /// Quality scores (FASTQ only; empty for FASTA).
    pub fn quality_scores(&self) -> &[u8] {
        &self.qual
    }

    /// Whether the record came from a FASTQ (`@`) header.
    pub fn is_fastq(&self) -> bool {
        self.is_fastq
    }
}

/// FAFQ sequential reader over any buffered input.
pub struct SeqReader<'a> {
    inner: Box<dyn BufRead + Send + 'a>,
    last: u8,
}

impl<'a> SeqReader<'a> {
    /// Opens a FAFQ reader from a path (supports stdin and gzip).
    pub fn new(infile: &str) -> anyhow::Result<SeqReader<'static>> {
        let inner: Box<dyn BufRead + Send> = crate::libs::io::reader(infile)
            .with_context(|| format!("Failed to open reader for {infile}"))?;
        Ok(SeqReader { inner, last: 0 })
    }

    /// Wraps an existing buffered reader.
    pub fn from_reader(inner: Box<dyn BufRead + Send + 'a>) -> Self {
        Self { inner, last: 0 }
    }

    /// Reads the next record into `rec`; returns `Ok(false)` at EOF.
    pub fn read_record(&mut self, rec: &mut SeqRecord) -> anyhow::Result<bool> {
        if self.last == 0 {
            match skip_to_header(&mut self.inner)? {
                Some(c) => self.last = c,
                None => return Ok(false),
            }
        }
        // Consume the header byte (`>`/`@`) at the front of the buffer: it
        // was left there either by `skip_to_header` or the previous
        // `read_seq` terminator.
        let buf = self.inner.fill_buf()?;
        if buf.first() == Some(&self.last) {
            self.inner.consume(1);
        }
        rec.clear();
        let term = read_definition(&mut self.inner, rec, self.last)?;
        self.last = if term == b'>' || term == b'@' {
            term
        } else {
            0
        };
        rec.is_fastq = term == b'+';
        if rec.is_fastq {
            let mut tmp = Vec::new();
            read_until_nl(&mut self.inner, &mut tmp)?;
            while rec.qual.len() < rec.seq.len() {
                if !read_until_nl(&mut self.inner, &mut rec.qual)? {
                    break; // EOF: truncated quality string
                }
            }
            if rec.qual.len() != rec.seq.len() {
                anyhow::bail!(
                    "FASTQ quality length {} != sequence length {}",
                    rec.qual.len(),
                    rec.seq.len()
                );
            }
        }
        Ok(true)
    }
}

/// Skips to the next `>` or `@`, returning the header byte (None at EOF).
fn skip_to_header<R: BufRead>(r: &mut R) -> std::io::Result<Option<u8>> {
    loop {
        let (found, consumed) = {
            let buf = r.fill_buf()?;
            if buf.is_empty() {
                return Ok(None);
            }
            match memchr2(b'>', b'@', buf) {
                Some(i) => (Some(buf[i]), i),
                None => (None, buf.len()),
            }
        };
        r.consume(consumed);
        if found.is_some() {
            // The header byte stays in the buffer; `read_record` consumes it.
            return Ok(found);
        }
    }
}

/// Reads name (to whitespace) and optional comment (to end of line), then the
/// sequence; returns the terminator (`>`, `@`, `+`) or 0 at EOF.
fn read_definition<R: BufRead>(r: &mut R, rec: &mut SeqRecord, header: u8) -> std::io::Result<u8> {
    let _ = header;
    let delim = read_until_ws(r, &mut rec.name)?;
    if delim != b'\n' && delim != b'\r' && delim != 0 {
        read_until_nl(r, &mut rec.comment)?;
    }
    read_seq(r, &mut rec.seq)
}

/// Appends bytes up to the next whitespace (` `, `\t`, `\n`, `\r`) and
/// returns the delimiter (0 at EOF); the delimiter is consumed.
fn read_until_ws<R: BufRead>(r: &mut R, out: &mut Vec<u8>) -> std::io::Result<u8> {
    loop {
        let (delim, consumed) = {
            let buf = r.fill_buf()?;
            if buf.is_empty() {
                return Ok(0);
            }
            match buf
                .iter()
                .position(|&b| b == b' ' || b == b'\t' || b == b'\n' || b == b'\r')
            {
                Some(i) => {
                    out.extend_from_slice(&buf[..i]);
                    (buf[i], i + 1)
                }
                None => {
                    out.extend_from_slice(buf);
                    (0, buf.len())
                }
            }
        };
        r.consume(consumed);
        if delim != 0 {
            return Ok(delim);
        }
    }
}

/// Appends one line (without `\n`, stripping a trailing `\r`); returns
/// whether a newline terminated the line (false at EOF).
fn read_until_nl<R: BufRead>(r: &mut R, out: &mut Vec<u8>) -> std::io::Result<bool> {
    loop {
        let (done, consumed) = {
            let buf = r.fill_buf()?;
            if buf.is_empty() {
                return Ok(false);
            }
            match memchr(b'\n', buf) {
                Some(i) => {
                    let line = buf[..i].strip_suffix(b"\r").unwrap_or(&buf[..i]);
                    out.extend_from_slice(line);
                    (true, i + 1)
                }
                None => {
                    out.extend_from_slice(buf);
                    (false, buf.len())
                }
            }
        };
        r.consume(consumed);
        if done {
            return Ok(true);
        }
    }
}

/// Appends sequence lines (newlines stripped) until a line starts with
/// `>`, `@` or `+`; returns the terminator (0 at EOF).
fn read_seq<R: BufRead>(r: &mut R, out: &mut Vec<u8>) -> std::io::Result<u8> {
    loop {
        let (term, consumed) = {
            let buf = r.fill_buf()?;
            if buf.is_empty() {
                return Ok(0);
            }
            if buf[0] == b'>' || buf[0] == b'@' || buf[0] == b'+' {
                return Ok(buf[0]);
            }
            match memchr(b'\n', buf) {
                Some(i) => {
                    let line = buf[..i].strip_suffix(b"\r").unwrap_or(&buf[..i]);
                    out.extend_from_slice(line);
                    (0, i + 1)
                }
                None => {
                    out.extend_from_slice(buf);
                    (0, buf.len())
                }
            }
        };
        r.consume(consumed);
        if term != 0 {
            return Ok(term);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_fasta as fasta;
    use noodles_fastq as fastq;
    use rand::{rngs::StdRng, Rng, SeedableRng};
    use std::io::Cursor;

    fn random_dna(rng: &mut StdRng, len: usize) -> Vec<u8> {
        let bases = b"ACGT";
        (0..len).map(|_| bases[rng.random_range(0..4)]).collect()
    }

    fn random_fasta(rng: &mut StdRng, n_seq: usize) -> Vec<u8> {
        let mut out = Vec::new();
        for s in 0..n_seq {
            let len = rng.random_range(1..500);
            let seq = random_dna(rng, len);
            out.extend_from_slice(format!(">seq{s} some comment\n").as_bytes());
            let mut col = 0usize;
            for &b in &seq {
                out.push(b);
                col += 1;
                if col == rng.random_range(40..100) {
                    out.push(b'\n');
                    col = 0;
                }
            }
            if col != 0 {
                out.push(b'\n');
            }
        }
        out
    }

    fn read_all_fa(data: Vec<u8>) -> Vec<(Vec<u8>, Vec<u8>)> {
        let mut r = SeqReader::from_reader(Box::new(Cursor::new(data)) as Box<dyn BufRead + Send>);
        let mut rec = SeqRecord::new();
        let mut out = Vec::new();
        while r.read_record(&mut rec).unwrap() {
            out.push((rec.name().to_vec(), rec.sequence().to_vec()));
        }
        out
    }

    fn read_all_fq(data: Vec<u8>) -> Vec<(Vec<u8>, Vec<u8>, Vec<u8>)> {
        let mut r = SeqReader::from_reader(Box::new(Cursor::new(data)) as Box<dyn BufRead + Send>);
        let mut rec = SeqRecord::new();
        let mut out = Vec::new();
        while r.read_record(&mut rec).unwrap() {
            out.push((
                rec.name().to_vec(),
                rec.sequence().to_vec(),
                rec.quality_scores().to_vec(),
            ));
        }
        out
    }

    #[test]
    fn fasta_matches_noodles_random() {
        let mut rng = StdRng::seed_from_u64(20260809);
        for _ in 0..10 {
            let fa = random_fasta(&mut rng, 30);
            let mut nr = fasta::io::Reader::new(&fa[..]);
            let mut expected = Vec::new();
            for rec in nr.records() {
                let rec = rec.unwrap();
                expected.push((rec.name().to_vec(), rec.sequence().as_ref().to_vec()));
            }
            assert_eq!(read_all_fa(fa), expected);
        }
    }

    #[test]
    fn fasta_single_line_and_crlf() {
        // Single-line records, CRLF line endings.
        let data = b">a\r\nACGTACGT\r\n>b\r\nTTTT\r\n>c\r\n";
        let got = read_all_fa(data.to_vec());
        assert_eq!(
            got,
            vec![
                (b"a".to_vec(), b"ACGTACGT".to_vec()),
                (b"b".to_vec(), b"TTTT".to_vec()),
                (b"c".to_vec(), b"".to_vec()),
            ]
        );
    }

    #[test]
    fn fasta_empty_lines_skipped() {
        let data = b">a\nACGT\n\nNN\n>b\n\n";
        let got = read_all_fa(data.to_vec());
        assert_eq!(
            got,
            vec![
                (b"a".to_vec(), b"ACGTNN".to_vec()),
                (b"b".to_vec(), b"".to_vec()),
            ]
        );
    }

    #[test]
    fn fastq_single_line_matches_noodles() {
        let data = b"@r1 desc\nACGTACGT\n+\nIIIIIIII\n@r2\nTTTT\n+\n!!!!\n";
        let mut nr = fastq::io::Reader::new(&data[..]);
        let mut expected = Vec::new();
        for rec in nr.records() {
            let rec = rec.unwrap();
            let seq: &[u8] = rec.sequence();
            let qual: &[u8] = rec.quality_scores();
            expected.push((rec.name().to_vec(), seq.to_vec(), qual.to_vec()));
        }
        assert_eq!(read_all_fq(data.to_vec()), expected);
    }

    #[test]
    fn fastq_multiline_supported() {
        // Multi-line seq/qual: noodles_fastq records() rejects this, but the
        // self-contained reader must accept it (kseq semantics).
        let data = b"@r1\nACGT\nACGT\n+\nIIII\nIIII\n";
        let got = read_all_fq(data.to_vec());
        assert_eq!(
            got,
            vec![(b"r1".to_vec(), b"ACGTACGT".to_vec(), b"IIIIIIII".to_vec(),)]
        );
    }

    #[test]
    fn fastq_quality_mismatch_errors() {
        let data = b"@r1\nACGT\n+\nIII\n";
        let mut r = SeqReader::from_reader(Box::new(&data[..]));
        let mut rec = SeqRecord::new();
        assert!(r.read_record(&mut rec).is_err());
    }

    #[test]
    fn mixed_fafq_auto_detect() {
        // One stream with FASTA then FASTQ records; the reader must
        // auto-detect each header.
        let data = b">fa1\nACGT\n@fq1\nTTTT\n+\n!!!!\n>fa2\nGGGG\n";
        let mut r = SeqReader::from_reader(Box::new(&data[..]));
        let mut rec = SeqRecord::new();
        assert!(r.read_record(&mut rec).unwrap());
        assert!(!rec.is_fastq());
        assert_eq!(rec.sequence(), b"ACGT");
        assert!(r.read_record(&mut rec).unwrap());
        assert!(rec.is_fastq());
        assert_eq!(rec.sequence(), b"TTTT");
        assert_eq!(rec.quality_scores(), b"!!!!");
        assert!(r.read_record(&mut rec).unwrap());
        assert!(!rec.is_fastq());
        assert_eq!(rec.sequence(), b"GGGG");
        assert!(!r.read_record(&mut rec).unwrap());
    }

    #[test]
    fn non_utf8_name_is_byte_clean() {
        // bstr semantics: the name is a byte string; non-UTF-8 bytes must not
        // fail at the read layer (consumers choose how to decode).
        let data = b">seq\xff\xfe desc\nACGT\n";
        let mut r = SeqReader::from_reader(Box::new(Cursor::new(data.to_vec())));
        let mut rec = SeqRecord::new();
        assert!(r.read_record(&mut rec).unwrap());
        assert_eq!(rec.name(), b"seq\xff\xfe");
        assert_eq!(rec.comment(), b"desc");
        assert_eq!(rec.sequence(), b"ACGT");
    }

    #[test]
    fn comment_and_name_parsing() {
        let data = b">id1 desc1 desc2\nACGT\n";
        let mut r = SeqReader::from_reader(Box::new(&data[..]));
        let mut rec = SeqRecord::new();
        r.read_record(&mut rec).unwrap();
        assert_eq!(rec.name(), b"id1");
        assert_eq!(rec.comment(), b"desc1 desc2");
        assert_eq!(rec.sequence(), b"ACGT");
    }
}
