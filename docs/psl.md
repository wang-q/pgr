# pgr psl

`pgr psl` provides tools for manipulating **PSL** alignment files (UCSC format).

## Subcommands

*   `chain`: Chain PSL alignments (connect alignment blocks).
*   `histo`: Collect alignment statistics for histograms.
*   `lift`: Lift PSL coordinates from fragment alignments to genomic coordinates.
*   `rc`: Reverse-complement alignments.
*   `stats`: Collect statistics from a PSL file (per-query or per-alignment).
*   `swap`: Swap target and query.
*   `to-chain`: Convert PSL to Chain format.
*   `to-range`: Extract alignment coordinates as ranges (.rg).

---

## chain

Connects alignment blocks in a PSL file into "chains" using dynamic programming. This is similar to UCSC `axtChain` but works directly on PSL files.

```bash
pgr psl chain [OPTIONS] <target.2bit> <query.2bit> <in.psl>
```

### Options

*   `--linear-gap <loose|medium>`: Linear gap cost presets (default: `loose`).
*   `--gap-open <int>`: Gap open cost (overrides linear-gap).
*   `--gap-extend <int>`: Gap extension cost (overrides linear-gap).
*   `--score-scheme <file|preset>`: Scoring matrix (LASTZ format) or preset (e.g., `hoxd55`).
*   `--min-score <float>`: Minimum chain score to output (default: 1000).
*   `-o, --outfile <file>`: Output filename (default: stdout).

### Examples

1.  **Chain PSL alignments with default settings**:
    ```bash
    pgr psl chain t.2bit q.2bit in.psl -o out.chain
    ```

2.  **Chain with affine gap costs**:
    ```bash
    pgr psl chain t.2bit q.2bit in.psl -o out.chain --gap-open 400 --gap-extend 30
    ```

---

## histo

Collects counts on PSL alignments for making histograms (e.g., alignments per query, coverage spread).

```bash
pgr psl histo [OPTIONS] --what <TYPE> <input>
```

### Options

*   `--what <alignsPerQuery|coverSpread|idSpread>`: Data to collect.
    *   `alignsPerQuery`: Number of alignments per query.
    *   `coverSpread`: Difference between max and min coverage for a query.
    *   `idSpread`: Difference between max and min identity for a query.
*   `--multi-only`: Omit queries with only one alignment.
*   `--non-zero`: Omit queries with zero values.
*   `-o, --output <file>`: Output filename (default: stdout).

### Examples

1.  **Count alignments per query**:
    ```bash
    pgr psl histo --what alignsPerQuery in.psl -o out.histo
    ```

---

## lift

Lifts PSL coordinates from fragment alignments (e.g., `chr1:100-200`) back to genomic coordinates.

```bash
pgr psl lift [OPTIONS] <infile>
```

### Options

*   `--q-sizes <file>`: Query sizes file (name, size).
*   `--t-sizes <file>`: Target sizes file (name, size).
*   `-o, --outfile <file>`: Output filename (default: stdout).

### Examples

1.  **Lift query coordinates**:
    ```bash
    pgr psl lift input.psl --q-sizes chrom.sizes > output.psl
    ```

---

## rc

Reverse-complements alignments in a PSL file.

```bash
pgr psl rc [OPTIONS] <input>
```

### Examples

1.  **Reverse complement PSL**:
    ```bash
    pgr psl rc in.psl -o out.psl
    ```

---

## stats

Collects statistics from a PSL file. Can report per-alignment, per-query, or overall statistics.

```bash
pgr psl stats [OPTIONS] <input>
```

### Options

*   `--query-stats`: Output per-query statistics.
*   `--overall-stats`: Output overall statistics.
*   `--queries <file>`: Tab-separated file with expected qNames and sizes.
*   `--tsv`: Write TSV header instead of autoSql style header.
*   `-o, --output <file>`: Output filename (default: stdout).

### Examples

1.  **Per-alignment statistics (default)**:
    ```bash
    pgr psl stats in.psl -o out.stats
    ```

2.  **Per-query statistics**:
    ```bash
    pgr psl stats --query-stats in.psl -o out.stats
    ```

---

## swap

Swaps target and query in a PSL file.

```bash
pgr psl swap [OPTIONS] <input>
```

### Options

*   `--no-rc`: Don't reverse-complement; just swap and make target strand explicit (if needed).
*   `-o, --output <file>`: Output filename (default: stdout).

### Examples

1.  **Swap target and query**:
    ```bash
    pgr psl swap in.psl -o out.psl
    ```

---

## to-chain

Converts PSL format to Chain format.

```bash
pgr psl to-chain [OPTIONS] <input>
```

### Options

*   `--fix-strand`: Fix `-` target strand by reverse complementing the record.
*   `-o, --output <file>`: Output filename (default: stdout).

### Examples

1.  **Convert PSL to Chain**:
    ```bash
    pgr psl to-chain in.psl -o out.chain
    ```

---

## to-range

Extracts alignment coordinates from PSL files as ranges (chr:start-end).

```bash
pgr psl to-range [OPTIONS] <infile>
```

### Options

*   `--target`: Extract target coordinates instead of query.
*   `-o, --outfile <file>`: Output filename (default: stdout).

### Examples

1.  **Extract query ranges**:
    ```bash
    pgr psl to-range input.psl > query.rg
    ```

2.  **Extract target ranges**:
    ```bash
    pgr psl to-range input.psl --target > target.rg
    ```
