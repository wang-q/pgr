# pgr lav

`pgr lav` provides tools for manipulating **LAV** (Local Alignment View) files and wrapping the **LASTZ** aligner.

## Subcommands

*   `lastz`: A wrapper for LASTZ alignment (Cactus style).
*   `to-psl`: Convert LAV files to PSL format.

---

## lastz

This command wraps `lastz` to perform alignments suitable for the Cactus RepeatMasking workflow or general pairwise alignment tasks.

It handles parallel execution for multiple target/query files, directory recursion, and provides predefined parameter sets (presets) for common species comparisons.

**Note**: `lastz` must be installed and available in your `PATH`.

```bash
pgr lav lastz [OPTIONS] <target> <query>
```

### Options

*   `--depth <int>`: Query depth threshold (default: 50). This sets `--querydepth=keep,nowarn:N` for `lastz`, which stops processing a query sequence once the coverage depth exceeds N (approx. 50x coverage).
*   `--self`: Perform self-alignment (align target against itself). In this mode, `target` and `query` should point to the same file or directory.
*   `-s, --preset <set01..set07>`: Use a predefined parameter set (see below).
*   `--show-preset`: Display the configuration (parameters & matrix) for the selected preset and exit.
*   `--lastz-args <string>`: Additional arguments passed directly to `lastz` (overrides preset settings).
*   `-o, --output <dir>`: Output directory (default: "lastz_out").
*   `-p, --parallel <int>`: Number of parallel threads (default: 4).

### Presets

*   `set01`: Hg17vsPanTro1 (Human vs Chimp)
*   `set02`: Hg19vsPanTro2 (Human vs Primate, more sensitive)
*   `set03`: Hg17vsMm5 (Human vs Mouse)
*   `set04`: Hg17vsRheMac2 (Human vs Macaque)
*   `set05`: Hg17vsBosTau2 (Human vs Cow)
*   `set06`: Hg17vsDanRer3 (Human vs Zebrafish)
*   `set07`: Hg17vsMonDom1 (Human vs Opossum)

### Examples

1.  **Align single files using a preset**:
    ```bash
    pgr lav lastz target.fa query.fa --preset set01 -o lastz_out
    ```

2.  **Align all FASTA files in directories**:
    ```bash
    pgr lav lastz target_dir/ query_dir/ --preset set03 -o lastz_out
    ```

3.  **Show parameters for a preset**:
    ```bash
    pgr lav lastz --preset set01 --show-preset
    ```

---

## to-psl

Converts BLASTZ/LASTZ LAV format files to PSL format.

```bash
pgr lav to-psl [OPTIONS] [input]
```

### Options

*   `[input]`: Input LAV file (default: stdin).
*   `-o, --output <file>`: Output PSL file (default: stdout).
*   `--target-strand <strand>`: Set the target strand (e.g., "+"). Default is no strand info.

### Examples

1.  **Convert LAV to PSL**:
    ```bash
    pgr lav to-psl in.lav -o out.psl
    ```
