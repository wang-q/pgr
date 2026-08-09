use criterion::{black_box, criterion_group, criterion_main, Criterion};
use noodles_fasta as fasta;
use std::io::BufRead;

/// Deterministic 50 MB FASTA: 50 sequences x 1 Mbp, 80 bp lines.
fn generate_fasta(n_seq: usize, seq_len: usize) -> Vec<u8> {
    let bases = b"ACGT";
    let mut x = 0x1234_5678_9abc_def0u64;
    let mut out = Vec::with_capacity(n_seq * (seq_len + seq_len / 80 + 20));
    for s in 0..n_seq {
        out.extend_from_slice(format!(">seq{s}\n").as_bytes());
        let mut col = 0usize;
        for _ in 0..seq_len {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            out.push(bases[(x >> 33) as usize & 3]);
            col += 1;
            if col == 80 {
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

fn bench_fa_read(c: &mut Criterion) {
    let fasta = generate_fasta(50, 1_000_000);
    let mut group = c.benchmark_group("fa_read_50mb");
    group.sample_size(10);

    group.bench_function("noodles_records", |b| {
        b.iter(|| {
            let mut reader = fasta::io::Reader::new(&fasta[..]);
            let mut total = 0usize;
            for rec in reader.records() {
                let rec = rec.unwrap();
                total += rec.sequence().len();
            }
            black_box(total)
        })
    });
    group.bench_function("noodles_sequence_only", |b| {
        b.iter(|| {
            // Lower-level: iterate definitions + read sequence bytes, no
            // Record construction.
            let mut cursor = &fasta[..];
            let mut total = 0usize;
            loop {
                let mut def = Vec::new();
                let n = read_until_nl(&mut cursor, &mut def);
                if n == 0 {
                    break;
                }
                let mut seq = Vec::new();
                read_sequence_until_def(&mut cursor, &mut seq);
                total += seq.len();
            }
            black_box(total)
        })
    });
    group.bench_function("naive_line_count", |b| {
        b.iter(|| {
            let mut cursor = &fasta[..];
            let mut total = 0usize;
            let mut line = Vec::new();
            loop {
                line.clear();
                let n = read_until_nl(&mut cursor, &mut line);
                if n == 0 {
                    break;
                }
                if line.first() != Some(&b'>') {
                    total += line.len();
                }
            }
            black_box(total)
        })
    });
    group.finish();
}

fn read_until_nl<R: BufRead>(r: &mut R, buf: &mut Vec<u8>) -> usize {
    let mut total = 0;
    loop {
        let (done, used) = {
            let available = match r.fill_buf() {
                Ok(b) => b,
                Err(_) => return total,
            };
            match available.iter().position(|&b| b == b'\n') {
                Some(i) => {
                    buf.extend_from_slice(&available[..i]);
                    (true, i + 1)
                }
                None => {
                    buf.extend_from_slice(available);
                    (false, available.len())
                }
            }
        };
        r.consume(used);
        total += used;
        if done || total == 0 {
            break;
        }
        if used == 0 {
            break;
        }
    }
    total
}

fn read_sequence_until_def<R: BufRead>(r: &mut R, buf: &mut Vec<u8>) {
    loop {
        let available = match r.fill_buf() {
            Ok(b) => b,
            Err(_) => return,
        };
        if available.is_empty() || available[0] == b'>' {
            return;
        }
        let i = match available.iter().position(|&b| b == b'\n') {
            Some(i) => i,
            None => available.len(),
        };
        buf.extend_from_slice(&available[..i]);
        r.consume(i + 1);
    }
}

criterion_group!(benches, bench_fa_read);
criterion_main!(benches);
