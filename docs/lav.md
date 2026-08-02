# pgr lav

`pgr lav` provides tools for manipulating **LAV** (Local Alignment View)
files. The LASTZ aligner wrapper lives under `pgr align`; see
[`align-lastz.md`](align-lastz.md).

## to-psl

Converts BLASTZ/LASTZ LAV format files to PSL format.

```bash
pgr lav to-psl [OPTIONS] [input]
```

### Options

*   `[input]`: Input LAV file (default: stdin).
*   `-o, --outfile <file>`: Output PSL file (default: stdout).
*   `--target-strand <strand>`: Set the target strand (e.g., "+"). Default is no strand info.
*   `--score-file <file>`: Output lav scores to side file (not yet implemented).

### Examples

1.  **Convert LAV to PSL**:
    ```bash
    pgr lav to-psl in.lav -o out.psl
    ```
