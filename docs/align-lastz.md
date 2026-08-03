# pgr align lastz

`pgr align lastz` wraps the external **LASTZ** aligner (Cactus style) to align
two genomes or directories of FASTA files:

```bash
pgr align lastz [OPTIONS] <target> <query>
```

It handles parallel execution for multiple target/query files, directory
recursion, and provides predefined parameter sets (presets) for common species
comparisons. The output is LAV format; convert it with `pgr lav to-psl`.

**Note**: `lastz` must be installed and available in your `PATH`.
**Note**: each input FASTA file must contain a single sequence (LAV output is
incompatible with lastz's `[multiple]` action) and must be plain text —
gzipped files are not read by lastz. Split multi-contig genomes with
`pgr fa split name` or pass directories of single-contig files.

## Options

*   `--query-depth <int>`: Query depth threshold (default: 50). This sets
    `--querydepth=keep,nowarn:N` for `lastz`, which stops processing a query
    sequence once the coverage depth exceeds N (approx. 50x coverage).
*   `--self`: Self-alignment (align the target against itself); the query
    input may be omitted or must be the same path as the target. Omitting the
    query implies `--self`.
*   `--preset <set01..set07>`: Use a predefined parameter set (see below).
*   `--show-preset`: Display the configuration (parameters & matrix) for the
    selected preset and exit.
*   `--lastz-args <string>`: Additional arguments passed directly to `lastz`
    (overrides preset settings).
*   `-o, --outdir <dir>`: Output directory (default: "lastz_out").
*   `-p, --parallel <int>`: Number of parallel threads (default: 4).

## Presets

*   `set01`: Hg17vsPanTro1 (Human vs Chimp)
*   `set02`: Hg19vsPanTro2 (Human vs Primate, more sensitive)
*   `set03`: Hg17vsMm5 (Human vs Mouse)
*   `set04`: Hg17vsRheMac2 (Human vs Macaque)
*   `set05`: Hg17vsBosTau2 (Human vs Cow)
*   `set06`: Hg17vsDanRer3 (Human vs Zebrafish)
*   `set07`: Hg17vsMonDom1 (Human vs Opossum)

## Examples

1.  **Align single files using a preset**:
    ```bash
    pgr align lastz target.fa query.fa --preset set01 -o lastz_out
    ```

2.  **Align all FASTA files in directories**:
    ```bash
    pgr align lastz target_dir/ query_dir/ --preset set03 -o lastz_out
    ```

3.  **Show parameters for a preset**:
    ```bash
    pgr align lastz --preset set01 --show-preset
    ```
4.  **Self-alignment with a single input**:
    ```bash
    pgr align lastz genome.fa --outdir self_out
    ```
