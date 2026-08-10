# CLI benchmark: pgr vs BBTools 39.38

Hyperfine 1.19.0 comparison of the pgr-native replacements against the
locally installed BBTools 39.38 (cbp package) on the Lambda test data
(40,000 reads, SRR5042715).

## Methodology

* `hyperfine --warmup 2 --runs 6` per pair; mean ± σ reported.
* pgr: parallel `target/release/pgr`; `fq trim-adapter` runs with
  `--threads 8` (matching BBTools); `fq clump`, `kmer hist`, and `fq norm`
  parallelize internally via rayon (logical CPU count, 32 on this host).
  BBTools: `threads=8`.
* Both sides write plain-text FASTQ (pgr does not gzip output); inputs are
  the same gzipped files (raw R1/R2, or the golden clumpify/trim/filter
  intermediates for the downstream steps).
* Parameters match the anchr trim pipeline (`k=31` clumpify, `k=23` trim,
  `k=27` filter, `k=31` khist; `ordered=t`, fixed seeds).
* Updated 2026-08-10 after pgr parallelization; the previous run was
  single-threaded.

## Results

| Step | BBTools 39.38 | pgr | Speedup |
|---|---:|---:|---:|
| clumpify / `fq clump` | 212.9 ± 5.0 ms | 88.1 ± 4.0 ms | 2.42x |
| bbduk trim / `fq trim-adapter` | 389.3 ± 47.9 ms | 80.6 ± 8.0 ms | 4.83x |
| bbduk filter / filter mode | 187.8 ± 6.4 ms | 74.4 ± 8.4 ms | 2.53x |
| kmercountexact / `kmer hist` | 2.554 ± 0.014 s | 682.5 ± 4.3 ms | 3.74x |
| repair / `fq split` | 174.8 ± 4.8 ms | 21.4 ± 0.7 ms | 8.16x |
| reformat sample / `fq sample` | 157.5 ± 3.0 ms | 24.0 ± 0.6 ms | 6.57x |
| bbnorm / `fq norm` | 838.7 ± 45.1 ms | 232.0 ± 6.5 ms | 3.62x |

Sum of step means: BBTools ≈ 4.51 s, pgr ≈ 1.20 s (~3.8x end to end, without
the JVM/script startup and intermediate gzip files of the real pipeline).

## Output equivalence

`clump`, `trim`, `filter`, `split`, `sample`, and the khist/peaks text are
byte-identical between the two sides on these inputs. `fq norm` keeps 39,846
reads vs bbnorm's 39,888 (exact canonical counts vs `bits=16` approximate
hash counts; see tests/bbtools/Lambda/README.md).

## Caveats

* `fq trim-adapter` is the main worker-pool parallelized command (`--threads`,
  default logical CPU count); on 50万-pair synthetic data it scaled from
  9.2 s (1 thread) to 1.4 s (8 threads), 6.6x, with ~15 MB peak memory.
* `fq clump`, `kmer hist`, and `fq norm` parallelize through rayon internally;
  the Lambda data is small, so their parallel gain is limited (clump even
  measured slightly slower than the earlier single-threaded run).
* `fq split` / `fq sample` are I/O-bound and remain single-threaded.
