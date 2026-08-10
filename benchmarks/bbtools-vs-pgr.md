# CLI benchmark: pgr vs BBTools 39.38

Hyperfine 1.19.0 comparison of the pgr-native replacements against the
locally installed BBTools 39.38 (cbp package) on the Lambda test data
(40,000 reads, SRR5042715).

## Methodology

* `hyperfine --warmup 2 --runs 6` per pair; mean ± σ reported.
* pgr: single-threaded `target/release/pgr`; BBTools: `threads=8`.
* Both sides write plain-text FASTQ (pgr does not gzip output); inputs are
  the same gzipped files (raw R1/R2, or the golden clumpify/trim/filter
  intermediates for the downstream steps).
* Parameters match the anchr trim pipeline (`k=31` clumpify, `k=23` trim,
  `k=27` filter, `k=31` khist; `ordered=t`, fixed seeds).

## Results

| Step | BBTools 39.38 | pgr | Speedup |
|---|---:|---:|---:|
| clumpify / `fq clump` | 225.4 ± 5.7 ms | 74.1 ± 1.6 ms | 3.04x |
| bbduk trim / `fq trim-adapter` | 408.5 ± 13.6 ms | 401.5 ± 3.3 ms | 1.02x |
| bbduk filter / filter mode | 196.5 ± 6.2 ms | 110.8 ± 1.8 ms | 1.77x |
| kmercountexact / `kmer hist` | 2.565 ± 0.015 s | 699.2 ± 7.6 ms | 3.67x |
| repair / `fq split` | 183.5 ± 3.4 ms | 21.7 ± 0.6 ms | 8.46x |
| reformat sample / `fq sample` | 198.4 ± 81.6 ms | 24.7 ± 0.5 ms | 8.02x |
| bbnorm / `fq norm` | 888.1 ± 20.0 ms | 239.0 ± 4.7 ms | 3.72x |

Sum of step means: BBTools ≈ 4.67 s, pgr ≈ 1.57 s (~3.0x end to end, without
the JVM/script startup and intermediate gzip files of the real pipeline).

## Output equivalence

`clump`, `trim`, `filter`, `split`, `sample`, and the khist/peaks text are
byte-identical between the two sides on these inputs. `fq norm` keeps 39,846
reads vs bbnorm's 39,888 (exact canonical counts vs `bits=16` approximate
hash counts; see tests/bbtools/Lambda/README.md).

## Caveats

* pgr is single-threaded; the trim step is dominated by the ported BBDuk
  algorithm itself (essentially tied at 1.02x), while JVM-dominated steps
  (clump, split, sample, khist) show the largest gains.
* The `reformat sample` run had outliers (161-365 ms); treat its speedup as
  approximate.
