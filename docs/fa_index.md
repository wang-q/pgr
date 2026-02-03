# Index

`samtools faidx` is optimized for massive reference genomes, relying on fixed line widths to
efficiently seek data without loading entire chromosomes. However, this requires strictly formatted
input, often causing errors with draft assemblies or "messy" files.

`pgr fa` redefines FASTA access with a focus on **robustness**:

*   **Format Agnostic**: Unlike `samtools`, `pgr` indexes *any* valid FASTA file regardless of line
    wrapping. `pgr fa` reliably extracts subsequences from draft assemblies and "messy" files
    without prior cleanup.
*   **Unified Architecture**: Uses a consistent indexing strategy for both plain text and
    BGZF-compressed data, abstracting away compression details for seamless access.
*   **Performance Optimization**: An internal LRU cache accelerates access for microbial genomes,
    contigs, and protein sequences, minimizing disk I/O overhead during intensive retrieval operations.
*   **Large Genome Support**: For mammalian-sized genomes, `pgr 2bit` ports the UCSC 2bit tools while
    maintaining a consistent command-line interface with `pgr fa`, ensuring a uniform experience
    across formats.

```bash
# gz
bgzip --keep --index tests/index/final.contigs.fa

pgr fa gz tests/index/final.contigs.fa -o tmp.gz

faToTwoBit tests/index/final.contigs.fa tests/index/final.contigs.2bit

# range
samtools faidx tests/index/final.contigs.fa
samtools faidx tests/index/final.contigs.fa \
    "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_170:1-20" "k81_158:70001-70020"
samtools faidx tests/index/final.contigs.fa -r tests/index/sample.rg

pgr fa range tests/index/final.contigs.fa.gz \
    "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_170(-):1-20" "k81_158:70001-70020"
pgr fa range tests/index/final.contigs.fa.gz -r tests/index/sample.rg

pgr 2bit range tests/index/final.contigs.2bit \
    "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_170(-):1-20" "k81_158:70001-70020"
pgr 2bit range tests/index/final.contigs.2bit -r tests/index/sample.rg

```

## Benchmarks

```shell
cbp install samtools
cbp install chainnet
cbp install hyperfine

```

### `bgzip`

```shell
hyperfine --warmup 5 --export-markdown gz.md.tmp \
    -n "bgzip" \
    'rm -f tests/index/final.contigs.fa.gz*;
     bgzip --keep --index tests/index/final.contigs.fa' \
    -n "bgzip --threads 2" \
    'rm -f tests/index/final.contigs.fa.gz*;
     bgzip --keep --index --threads 2 tests/index/final.contigs.fa' \
    -n "pgr fa gz" \
    'rm -f tests/index/final.contigs.fa.gz*;
     pgr fa gz tests/index/final.contigs.fa' \
    -n "pgr fa gz -p 2" \
    'rm -f tests/index/final.contigs.fa.gz*;
     pgr fa gz -p 2 tests/index/final.contigs.fa' \
    -n "faToTwoBit" \
    'rm -f tests/index/final.contigs.2bit;
     faToTwoBit tests/index/final.contigs.fa tests/index/final.contigs.2bit' \
    -n "pgr fa to-2bit" \
    'rm -f tests/index/final.contigs.2bit;
     pgr fa to-2bit tests/index/final.contigs.fa -o tests/index/final.contigs.2bit'

cat gz.md.tmp

```

| Command             |  Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:--------------------|-----------:|---------:|---------:|------------:|
| `bgzip`             | 71.6 ± 1.5 |     69.2 |     74.6 | 2.35 ± 0.09 |
| `bgzip --threads 2` | 52.2 ± 0.9 |     50.6 |     54.4 | 1.71 ± 0.06 |
| `pgr fa gz`         | 42.2 ± 0.9 |     40.8 |     44.6 | 1.38 ± 0.05 |
| `pgr fa gz -p 2`    | 32.7 ± 1.3 |     30.7 |     40.5 | 1.07 ± 0.06 |
| `faToTwoBit`        | 30.5 ± 1.0 |     29.0 |     33.7 |        1.00 |
| `pgr fa to-2bit`    | 31.1 ± 2.4 |     28.1 |     38.3 | 1.02 ± 0.08 |

### Create `.loc` and `.fai`

```shell
hyperfine --warmup 5 --export-markdown faidx.md.tmp \
    -n "samtools faidx .fa" \
    'rm -f tests/index/final.contigs.fa.fai;
     samtools faidx tests/index/final.contigs.fa' \
    -n "samtools faidx .fa.gz" \
    'rm -f tests/index/final.contigs.fa.gz.fai;
     samtools faidx tests/index/final.contigs.fa.gz' \
    -n "pgr fa range .fa" \
    'rm -f tests/index/final.contigs.fa.loc;
     pgr fa range tests/index/final.contigs.fa' \
    -n "pgr fa range .fa.gz" \
    'rm -f tests/index/final.contigs.fa.gz.loc;
     pgr fa range tests/index/final.contigs.fa.gz'

cat faidx.md.tmp

```

| Command                 |  Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:------------------------|-----------:|---------:|---------:|------------:|
| `samtools faidx .fa`    | 17.7 ± 1.0 |     16.3 |     21.7 |        1.00 |
| `samtools faidx .fa.gz` | 21.2 ± 0.9 |     19.6 |     23.2 | 1.20 ± 0.09 |
| `pgr fa range .fa`      | 20.8 ± 0.9 |     19.4 |     24.9 | 1.17 ± 0.08 |
| `pgr fa range .fa.gz`   | 19.1 ± 0.8 |     17.7 |     22.0 | 1.08 ± 0.08 |

### `pgr fa range`

```shell
hyperfine --warmup 5 --export-markdown rg.md.tmp \
    -n "samtools faidx .fa" \
    'samtools faidx tests/index/final.contigs.fa \
        "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_158:70001-70020"' \
    -n "samtools faidx .fa.gz" \
    'samtools faidx tests/index/final.contigs.fa.gz \
        "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_158:70001-70020"' \
    -n "pgr fa range .fa" \
    'pgr fa range tests/index/final.contigs.fa \
        "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_158:70001-70020"' \
    -n "pgr fa range .fa.gz" \
    'pgr fa range tests/index/final.contigs.fa.gz \
        "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_158:70001-70020"' \
    -n "pgr 2bit range" \
    'pgr 2bit range tests/index/final.contigs.2bit \
        "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_158:70001-70020"'

cat rg.md.tmp

```

| Command                 |  Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:------------------------|-----------:|---------:|---------:|------------:|
| `samtools faidx .fa`    |  6.1 ± 0.4 |      4.8 |      8.5 | 1.07 ± 0.11 |
| `samtools faidx .fa.gz` |  7.7 ± 0.7 |      7.0 |     18.1 | 1.36 ± 0.16 |
| `pgr fa range .fa`      |  8.3 ± 0.6 |      7.3 |     10.4 | 1.45 ± 0.15 |
| `pgr fa range .fa.gz`   | 11.2 ± 0.7 |      9.8 |     13.0 | 1.96 ± 0.19 |
| `pgr 2bit range`        |  5.7 ± 0.4 |      5.0 |      7.2 |        1.00 |

### `pgr fa range -r`

```shell
hyperfine --warmup 5 --export-markdown rg.md.tmp \
    -n "samtools faidx .fa" \
    'samtools faidx tests/index/final.contigs.fa -r tests/index/sample.rg > /dev/null' \
    -n "samtools faidx .fa.gz" \
    'samtools faidx tests/index/final.contigs.fa.gz -r tests/index/sample.rg > /dev/null' \
    -n "pgr fa range .fa" \
    'pgr fa range tests/index/final.contigs.fa -r tests/index/sample.rg > /dev/null' \
    -n "pgr fa range .fa.gz" \
    'pgr fa range tests/index/final.contigs.fa.gz -r tests/index/sample.rg > /dev/null' \
    -n "pgr 2bit range" \
    'pgr 2bit range tests/index/final.contigs.2bit -r tests/index/sample.rg > /dev/null'

cat rg.md.tmp

```

| Command                   |  Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:--------------------------|-----------:|---------:|---------:|------------:|
| `samtools faidx .fa`      |  8.3 ± 0.7 |      7.1 |     12.1 | 1.05 ± 0.11 |
| `samtools faidx .fa.gz`   |  9.8 ± 0.5 |      9.0 |     11.9 | 1.25 ± 0.10 |
| `pgr fa range .fa`        | 10.0 ± 0.5 |      9.2 |     13.3 | 1.27 ± 0.10 |
| `pgr fa range .fa.gz`     | 12.4 ± 0.5 |     11.5 |     15.5 | 1.58 ± 0.12 |
| `pgr 2bit range`          |  7.9 ± 0.5 |      7.1 |     10.6 |        1.00 |
