use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};
use pgr::libs::paf::cigar::{classify_alignment, scan_cigar_ops, scan_cs};
use rand::Rng;

// 95% match / 3% query gaps / 2% mismatches, same profile as the 40 M-column
// `maf to-paf` end-to-end measurement (bench-profile-hotspots.md).
fn random_aln(rng: &mut impl Rng, n: usize) -> (Vec<u8>, Vec<u8>) {
    let bases = *b"ACGT";
    let r: Vec<u8> = (0..n).map(|_| bases[rng.random_range(0..4)]).collect();
    let mut q: Vec<u8> = Vec::with_capacity(n);
    for &b in &r {
        let x = rng.random_range(0..100);
        q.push(if x < 3 {
            b'-'
        } else if x < 5 {
            bases[rng.random_range(0..4)]
        } else {
            b
        });
    }
    (r, q)
}

// Pre-optimization two-pass implementation (per-column classification), kept
// as the scalar reference for the function-level speedup measurement.
fn old_cigar_ops(r: &[u8], q: &[u8]) -> usize {
    let mut ops: Vec<(u32, char)> = Vec::new();
    for (&rc, &qc) in r.iter().zip(q.iter()) {
        let op_char = match (rc, qc) {
            (b'-', b'-') => continue,
            (b'-', _) => 'I',
            (_, b'-') => 'D',
            _ if rc.eq_ignore_ascii_case(&qc) => '=',
            _ => 'X',
        };
        match ops.last_mut() {
            Some(last) if last.1 == op_char => last.0 += 1,
            _ => ops.push((1, op_char)),
        }
    }
    ops.len()
}

fn old_cs(r: &[u8], q: &[u8]) -> String {
    let mut cs = String::new();
    let mut run = 0usize;
    for (&rc, &qc) in r.iter().zip(q.iter()) {
        match (rc, qc) {
            (b'-', b'-') => continue,
            (b'-', qq) => {
                if run > 0 {
                    cs.push(':');
                    cs.push_str(&run.to_string());
                    run = 0;
                }
                cs.push('+');
                cs.push(qq.to_ascii_uppercase() as char);
            }
            (rr, b'-') => {
                if run > 0 {
                    cs.push(':');
                    cs.push_str(&run.to_string());
                    run = 0;
                }
                cs.push('-');
                cs.push(rr.to_ascii_uppercase() as char);
            }
            (rr, qq) if rr.eq_ignore_ascii_case(&qq) => run += 1,
            (rr, qq) => {
                if run > 0 {
                    cs.push(':');
                    cs.push_str(&run.to_string());
                    run = 0;
                }
                cs.push('*');
                cs.push(rr.to_ascii_uppercase() as char);
                cs.push(qq.to_ascii_uppercase() as char);
            }
        }
    }
    if run > 0 {
        cs.push(':');
        cs.push_str(&run.to_string());
    }
    cs
}

fn bench_cigar(c: &mut Criterion) {
    let mut rng = rand::rng();
    for (name, n) in [("10m", 10_000_000usize), ("40m", 40_000_000usize)] {
        let (r, q) = random_aln(&mut rng, n);
        let mut group = c.benchmark_group(format!("cigar_{name}"));
        group.bench_function("old_two_pass", |b| {
            b.iter_batched(
                || (black_box(r.clone()), black_box(q.clone())),
                |(r, q)| {
                    let ops = old_cigar_ops(&r, &q);
                    let cs = old_cs(&r, &q);
                    black_box((ops, cs.len()))
                },
                BatchSize::LargeInput,
            )
        });
        group.bench_function("new_classify_scan", |b| {
            b.iter_batched(
                || (black_box(r.clone()), black_box(q.clone())),
                |(r, q)| {
                    let mask = classify_alignment(&r, &q).unwrap();
                    let ops = scan_cigar_ops(&mask);
                    let cs = scan_cs(&mask, &r, &q);
                    black_box((ops.len(), cs.len()))
                },
                BatchSize::LargeInput,
            )
        });
        group.bench_function("new_classify_only", |b| {
            b.iter_batched(
                || (black_box(r.clone()), black_box(q.clone())),
                |(r, q)| {
                    let mask = classify_alignment(&r, &q).unwrap();
                    black_box(mask.m.len())
                },
                BatchSize::LargeInput,
            )
        });
        group.finish();
    }
}

criterion_group!(benches, bench_cigar);
criterion_main!(benches);
