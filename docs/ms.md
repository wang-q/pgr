# pgr ms

`pgr ms` provides tools for working with **Hudson's ms** simulator output.

## Subcommands

*   `to-dna`: Convert ms output haplotypes (0/1) to DNA sequences (FASTA).

---

## to-dna

Converts the 0/1 haplotype output from `ms` simulations into actual DNA sequences (FASTA format) by simulating an ancestral sequence and applying mutations.

```bash
pgr ms to-dna [OPTIONS] [files]...
```

### Options

*   `-g, --gc <float>`: GC content ratio for the ancestral sequence (0.0 to 1.0, default: 0.5).
*   `-s, --seed <int>`: Random seed. If omitted, uses system time and PID.
*   `--no-perturb`: Disable position micro-perturbation. By default, `ms` positions (0..1) are mapped to integer sites, and `pgr` slightly perturbs them to avoid collisions.
*   `-v, --verbose`: Print runtime information (paths, inputs, seed).
*   `--doc`: Print full documentation (this help).
*   `-o, --outfile <file>`: Output filename (default: stdout).

### Input/Output

*   **Input**: `ms` output files. Reads from stdin if no files are provided.
*   **Output**: FASTA format with single-line sequences.
    *   Headers format: `>[Lx_][Px_]Sx`
        *   `Lx`: Batch/Replicate index (if multiple replicates).
        *   `Px`: Population index (if multiple populations).
        *   `Sx`: Sample index.

### Examples

1.  **Pipe `ms` output directly**:
    ```bash
    ms 10 1 -t 5 -r 0 1000 | pgr ms to-dna --gc 0.5 > out.fa
    ```

2.  **Convert an existing ms file**:
    ```bash
    pgr ms to-dna input.ms -o out.fa --seed 12345
    ```

3.  **Disable position perturbation**:
    ```bash
    pgr ms to-dna input.ms --no-perturb
    ```
