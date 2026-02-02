# Index

`samtools faidx` is optimized for massive reference genomes, relying on fixed line widths to
efficiently seek data without loading entire chromosomes. However, this requires strictly formatted
input, often causing errors with draft assemblies or "messy" files.

`pgr fa` redefines FASTA access with a focus on **robustness**:

*   **Format Agnostic**: Unlike `samtools`, `pgr` indexes *any* valid FASTA file regardless of line
    wrapping. `pgr fa` reliably extracts subsequences from draft assemblies and "messy" files
    without prior cleanup.
*   **Unified Architecture**: Uses a consistent indexing strategy for both plain text and
    BGZF-compressed data. The `range` command works identically on both, abstracting away
    compression details.
*   **Performance Optimization**: An internal LRU cache accelerates access for microbial genomes and
    contigs, minimizing disk I/O overhead during intensive retrieval operations.
*   **Large Genome Support**: For mammalian-sized genomes, `pgr 2bit` ports the UCSC 2bit tools while
    maintaining a consistent command-line interface with `pgr fa`, ensuring a uniform experience
    across formats.

```bash
# gz
bgzip --keep --index tests/index/final.contigs.fa

pgr fa gz tests/index/final.contigs.fa -o tmp.gz

# range
samtools faidx tests/index/final.contigs.fa
samtools faidx tests/index/final.contigs.fa \
    "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_170:1-20" "k81_158:70001-70020"
samtools faidx tests/index/final.contigs.fa -r tests/index/sample.rg

hnsm range tests/index/final.contigs.fa.gz \
    "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_170(-):1-20" "k81_158:70001-70020"
hnsm range tests/index/final.contigs.fa.gz -r tests/index/sample.rg

```

## Benchmarks

```shell
cbp install samtools
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
     pgr fa gz -p 2 tests/index/final.contigs.fa'

cat gz.md.tmp

```

| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `bgzip` | 71.2 ± 1.7 | 68.9 | 76.3 | 2.12 ± 0.10 |
| `bgzip --threads 2` | 51.5 ± 1.0 | 50.0 | 54.8 | 1.53 ± 0.07 |
| `pgr fa gz` | 42.3 ± 1.2 | 40.1 | 45.0 | 1.26 ± 0.06 |
| `pgr fa gz -p 2` | 33.6 ± 1.3 | 32.0 | 40.6 | 1.00 |

### `.loc` and `.fai`

```shell
hyperfine --warmup 5 --export-markdown faidx.md.tmp \
    -n "samtools faidx .fa" \
    'rm -f tests/index/final.contigs.fa.fai;
     samtools faidx tests/index/final.contigs.fa' \
    -n "samtools faidx .fa.gz" \
    'rm -f tests/index/final.contigs.fa.gz.fai;
     samtools faidx tests/index/final.contigs.fa.gz' \
    -n "hnsm range .fa" \
    'rm -f tests/index/final.contigs.fa.loc;
     hnsm range tests/index/final.contigs.fa' \
    -n "hnsm range .fa.gz" \
    'rm -f tests/index/final.contigs.fa.gz.loc;
     hnsm range tests/index/final.contigs.fa.gz'

cat faidx.md.tmp

```

| Command                 |  Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:------------------------|-----------:|---------:|---------:|------------:|
| `samtools faidx .fa`    | 22.9 ± 1.5 |     20.1 |     27.8 |        1.00 |
| `samtools faidx .fa.gz` | 30.3 ± 1.5 |     27.6 |     36.1 | 1.32 ± 0.11 |
| `hnsm range .fa`        | 25.7 ± 1.4 |     23.2 |     33.2 | 1.12 ± 0.10 |
| `hnsm range .fa.gz`     | 25.2 ± 1.6 |     23.0 |     35.9 | 1.10 ± 0.10 |

### `hnsm range`

```shell
hyperfine --warmup 5 --export-markdown rg.md.tmp \
    -n "samtools faidx .fa" \
    'samtools faidx tests/index/final.contigs.fa \
        "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_158:70001-70020"' \
    -n "samtools faidx .fa.gz" \
    'samtools faidx tests/index/final.contigs.fa.gz \
        "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_158:70001-70020"' \
    -n "hnsm range .fa" \
    'hnsm range tests/index/final.contigs.fa \
        "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_158:70001-70020"' \
    -n "hnsm range .fa.gz" \
    'hnsm range tests/index/final.contigs.fa.gz \
        "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_158:70001-70020"'

cat rg.md.tmp

```

| Command                 |  Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:------------------------|-----------:|---------:|---------:|------------:|
| `samtools faidx .fa`    | 12.1 ± 0.8 |     10.3 |     14.9 | 1.21 ± 0.11 |
| `samtools faidx .fa.gz` | 14.3 ± 0.8 |     12.6 |     17.1 | 1.42 ± 0.12 |
| `hnsm range .fa`        | 10.0 ± 0.7 |      8.7 |     12.8 |        1.00 |
| `hnsm range .fa.gz`     | 14.1 ± 0.7 |     12.0 |     15.8 | 1.40 ± 0.12 |

### `hnsm range -r`

```shell
hyperfine --warmup 5 --export-markdown rg.md.tmp \
    -n "samtools faidx .fa" \
    'samtools faidx tests/index/final.contigs.fa -r tests/index/sample.rg > /dev/null' \
    -n "samtools faidx .fa.gz" \
    'samtools faidx tests/index/final.contigs.fa.gz -r tests/index/sample.rg > /dev/null' \
    -n "hnsm range .fa" \
    'hnsm range tests/index/final.contigs.fa -r tests/index/sample.rg > /dev/null' \
    -n "hnsm range .fa.gz" \
    'hnsm range tests/index/final.contigs.fa.gz -r tests/index/sample.rg > /dev/null'

cat rg.md.tmp

```

| Command                 |  Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:------------------------|-----------:|---------:|---------:|------------:|
| `samtools faidx .fa`    | 13.7 ± 0.7 |     11.6 |     16.8 | 1.16 ± 0.09 |
| `samtools faidx .fa.gz` | 16.2 ± 0.9 |     14.4 |     19.1 | 1.38 ± 0.11 |
| `hnsm range .fa`        | 11.8 ± 0.7 |     10.3 |     13.9 |        1.00 |
| `hnsm range .fa.gz`     | 15.7 ± 0.8 |     14.2 |     18.8 | 1.33 ± 0.10 |
