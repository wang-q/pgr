use criterion::{black_box, criterion_group, criterion_main, Criterion};
use memchr::memchr;
use noodles_fasta as fasta;
use noodles_fastq as fastq;
use pgr::libs::fmt::seq::{SeqReader, SeqRecord};
use std::io::BufRead;

/// Deterministic 50 MB FASTA (50 x 1 Mbp, 80 bp lines) and a matching FASTQ
/// (same sequences, constant Phred+33 'I' qualities).
fn generate_fafq(n_seq: usize, seq_len: usize) -> (Vec<u8>, Vec<u8>) {
    let bases = b"ACGT";
    let mut x = 0x1234_5678_9abc_def0u64;
    let mut fa = Vec::with_capacity(n_seq * (seq_len + seq_len / 80 + 20));
    let mut fq = Vec::with_capacity(n_seq * (2 * seq_len + seq_len / 80 * 2 + 40));
    for s in 0..n_seq {
        let mut seq = Vec::with_capacity(seq_len);
        fa.extend_from_slice(format!(">seq{s}\n").as_bytes());
        let mut col = 0usize;
        for _ in 0..seq_len {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let b = bases[(x >> 33) as usize & 3];
            fa.push(b);
            seq.push(b);
            col += 1;
            if col == 80 {
                fa.push(b'\n');
                col = 0;
            }
        }
        if col != 0 {
            fa.push(b'\n');
        }
        // FASTQ: single-line seq + qual (noodles_fastq records() reads only
        // one line per field; multi-line FASTQ is not supported there).
        fq.extend_from_slice(format!("@seq{s}\n").as_bytes());
        fq.extend_from_slice(&seq);
        fq.push(b'\n');
        fq.extend_from_slice(b"+\n");
        fq.extend(std::iter::repeat_n(b'I', seq_len));
        fq.push(b'\n');
    }
    (fa, fq)
}

/// kseq-style record with reused buffers (capacity kept across reads).
#[derive(Default)]
struct KseqRec {
    name: Vec<u8>,
    comment: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    last: u8,
}

fn read_byte<R: BufRead>(r: &mut R) -> std::io::Result<Option<u8>> {
    let buf = r.fill_buf()?;
    if buf.is_empty() {
        return Ok(None);
    }
    let b = buf[0];
    r.consume(1);
    Ok(Some(b))
}

fn read_until_nl<R: BufRead>(r: &mut R, out: &mut Vec<u8>) -> std::io::Result<()> {
    loop {
        let (done, consumed) = {
            let buf = r.fill_buf()?;
            if buf.is_empty() {
                return Ok(());
            }
            match memchr(b'\n', buf) {
                Some(i) => {
                    out.extend_from_slice(&buf[..i]);
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
            return Ok(());
        }
    }
}

fn read_until_ws<R: BufRead>(r: &mut R, out: &mut Vec<u8>) -> std::io::Result<()> {
    loop {
        let (done, consumed) = {
            let buf = r.fill_buf()?;
            if buf.is_empty() {
                return Ok(());
            }
            match buf
                .iter()
                .position(|&b| b == b' ' || b == b'\t' || b == b'\n')
            {
                Some(i) => {
                    out.extend_from_slice(&buf[..i]);
                    (true, i)
                }
                None => {
                    out.extend_from_slice(buf);
                    (false, buf.len())
                }
            }
        };
        r.consume(consumed);
        if done {
            return Ok(());
        }
    }
}

/// Reads sequence lines until `>`, `@` or `+` at line start; returns the
/// terminator (0 at EOF). Newlines are not appended to `out`.
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
                    out.extend_from_slice(&buf[..i]);
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

/// kseq_read equivalent: one FAFQ record, buffers reused.
fn kseq_like_read<R: BufRead>(r: &mut R, rec: &mut KseqRec) -> std::io::Result<bool> {
    if rec.last == 0 {
        let mut c;
        loop {
            match read_byte(r)? {
                Some(b) => c = b,
                None => return Ok(false),
            }
            if c == b'>' || c == b'@' {
                break;
            }
        }
        rec.last = c;
    }
    rec.name.clear();
    rec.comment.clear();
    rec.seq.clear();
    rec.qual.clear();
    read_until_ws(r, &mut rec.name)?;
    read_until_nl(r, &mut rec.comment)?;
    let term = read_seq(r, &mut rec.seq)?;
    rec.last = if term == b'>' || term == b'@' {
        term
    } else {
        0
    };
    if term == b'+' {
        let mut tmp = Vec::new();
        read_until_nl(r, &mut tmp)?;
        while rec.qual.len() < rec.seq.len() {
            read_until_nl(r, &mut rec.qual)?;
        }
    }
    Ok(true)
}

fn bench_seq_read(c: &mut Criterion) {
    let (fa, fq) = generate_fafq(50, 1_000_000);
    let mut group = c.benchmark_group("seq_read_50mb");
    group.sample_size(10);

    group.bench_function("fasta_noodles_records", |b| {
        b.iter(|| {
            let mut reader = fasta::io::Reader::new(&fa[..]);
            let mut total = 0usize;
            for rec in reader.records() {
                total += rec.unwrap().sequence().len();
            }
            black_box(total)
        })
    });
    group.bench_function("fasta_kseq_like", |b| {
        b.iter(|| {
            let mut cursor = &fa[..];
            let mut rec = KseqRec::default();
            let mut total = 0usize;
            while kseq_like_read(&mut cursor, &mut rec).unwrap() {
                total += rec.seq.len();
            }
            black_box(total)
        })
    });
    group.bench_function("fastq_noodles_records", |b| {
        b.iter(|| {
            let mut reader = fastq::io::Reader::new(&fq[..]);
            let mut total = 0usize;
            for rec in reader.records() {
                total += rec.unwrap().sequence().len();
            }
            black_box(total)
        })
    });
    group.bench_function("fastq_kseq_like", |b| {
        b.iter(|| {
            let mut cursor = &fq[..];
            let mut rec = KseqRec::default();
            let mut total = 0usize;
            while kseq_like_read(&mut cursor, &mut rec).unwrap() {
                total += rec.seq.len();
            }
            black_box(total)
        })
    });
    group.bench_function("fasta_seq_reader_native", |b| {
        b.iter(|| {
            let mut r = SeqReader::from_reader(Box::new(std::io::Cursor::new(&fa[..])));
            let mut rec = SeqRecord::new();
            let mut total = 0usize;
            while r.read_record(&mut rec).unwrap() {
                total += rec.sequence().len();
            }
            black_box(total)
        })
    });
    group.bench_function("fastq_seq_reader_native", |b| {
        b.iter(|| {
            let mut r = SeqReader::from_reader(Box::new(std::io::Cursor::new(&fq[..])));
            let mut rec = SeqRecord::new();
            let mut total = 0usize;
            while r.read_record(&mut rec).unwrap() {
                total += rec.sequence().len();
            }
            black_box(total)
        })
    });
    group.finish();
}

criterion_group!(benches, bench_seq_read);
criterion_main!(benches);
