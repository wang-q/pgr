//! Benchmark: `Range` string parsing — regex decoder vs hand-written scanner.
//!
//! The corpus mixes the formats pgr actually parses (plain `chr:start-end`,
//! strand `chr(+):...`, species prefix `S288c.I(-):...`, single coordinates,
//! slash/underscore contigs) plus the regex fallback cases.

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pgr::libs::ds::Range;

/// The original regex decoder, kept here as the benchmark baseline (the
/// production `Range::from_str` is the hand-written scanner; the regex is
/// documented in `src/libs/ds/range.rs` and preserved as a test oracle).
fn regex_from_str(range: &str) -> Range {
    static RE: std::sync::LazyLock<regex::Regex> = std::sync::LazyLock::new(|| {
        regex::Regex::new(
            r"(?xi)
            (?:(?P<name>[\w_]+)\.)?
            (?P<chr>[\w/-]+)
            (?:\((?P<strand>.+)\))?
            [:]                    # spacer
            (?P<start>\d+)
            [_\-]?                 # spacer
            (?P<end>\d+)?
            ",
        )
        .expect("valid range regex")
    });
    let mut new = Range::new();
    let caps = match RE.captures(range) {
        Some(x) => x,
        None => {
            new.chr = range.split(' ').next().unwrap().to_string();
            return new;
        }
    };
    for name in RE.capture_names().flatten() {
        if let Some(m) = caps.name(name) {
            match name {
                "name" => new.name = m.as_str().to_string(),
                "chr" => new.chr = m.as_str().to_string(),
                "strand" => new.strand = m.as_str().to_string(),
                "start" => new.start = m.as_str().parse::<i32>().unwrap(),
                "end" => new.end = m.as_str().parse::<i32>().unwrap(),
                _ => {}
            }
        }
    }
    if new.start != 0 && new.end == 0 {
        new.end = new.start;
    }
    new
}

const CORPUS: &[&str] = &[
    "chr1:1-100",
    "chrM(+):1-16571",
    "NC_000913:100-200",
    "S288c.I(-):27070-29557",
    "S288c.II(+):1-813184",
    "1:1-23",
    "I:100",
    "I(+):100-200",
    "I(-):100-200",
    "infile_0/1/0_514:19-25",
    "a.b.c:1-2",
    "foo I:1-100",
    "S288c The baker's yeast",
    "1:-100",
    "chr1:1",
    "chr1:1_2",
    "chrX:1000000-2000000",
];

fn parse_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("range parse");
    group.bench_function("regex (baseline)", |b| {
        b.iter(|| {
            let mut acc = 0u64;
            for &s in CORPUS {
                let r = regex_from_str(black_box(s));
                acc = acc.wrapping_add(r.start as u64).wrapping_add(r.end as u64);
            }
            black_box(acc);
        })
    });
    group.bench_function("manual from_str", |b| {
        b.iter(|| {
            let mut acc = 0u64;
            for &s in CORPUS {
                let r = Range::from_str(black_box(s));
                acc = acc.wrapping_add(r.start as u64).wrapping_add(r.end as u64);
            }
            black_box(acc);
        })
    });
    group.finish();
}

criterion_group!(benches, parse_benchmark);
criterion_main!(benches);
