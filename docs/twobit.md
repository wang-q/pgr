# pgr 2bit

`pgr 2bit` provides tools for manipulating **2bit** files. 2bit is a binary format for storing genomic sequences efficiently (2 bits per base). It supports random access but does not support streaming (stdin) or gzip compression directly.

## Subcommands

The subcommands are organized into the following categories:

*   **Info**: Extract information or statistics.
    *   `masked`: Identify masked regions (soft or hard).
    *   `size`: Count total bases or sequence lengths.
*   **Subset**: Extract specific parts of the data.
    *   `range`: Extract sequence regions by coordinates.
    *   `some`: Extract full sequences based on a list of names.
*   **Transform**: Convert formats.
    *   `to-fa`: Convert 2bit to FASTA format.

---

## Info Commands

### masked

Identifies masked regions in one or more 2bit files.

```bash
pgr 2bit masked [OPTIONS] <infiles>...
```

*   `--gap`: Only identify hard-masked regions (N/n gaps).
*   `-o, --outfile <file>`: Output filename (default: stdout).
*   Output is one region per line in the form `seq_name:start-end` (1-based, inclusive). Single-base regions are emitted as `seq_name:position`.

### size

Retrieves sequence sizes (lengths) from one or more 2bit files.

```bash
pgr 2bit size [OPTIONS] <infiles>...
```

*   `--no-ns`: Output size excluding Ns (only A, C, G, T counts).
*   `-o, --outfile <file>`: Output filename (default: stdout).
*   Output is one `name\tsize` line per sequence, in the order sequences appear in each input file.

---

## Subset Commands

### range

Extracts sequence regions from 2bit files using genomic coordinates.

```bash
pgr 2bit range [OPTIONS] <infile> [ranges]...
```

*   `[ranges]`: List of ranges in format `seq_name(strand):start-end` or just `seq_name`.
    *   `seq_name` alone (without a colon) extracts the full sequence.
    *   `start-end` are 1-based, inclusive.
    *   `strand` is optional (`+` or `-`).
    *   Ranges extending past the end of a sequence are truncated with a warning.
*   `-r, --rgfile <file>`: File containing ranges, one per line.
*   `-o, --outfile <file>`: Output filename (default: stdout).
*   Output is FASTA; the header is the original range string (e.g. `>seq1:2-5` or `>seq1`).

### some

Extracts full sequences from a 2bit file based on a list of sequence names.

```bash
pgr 2bit some [OPTIONS] <infile> <list.txt>
```

*   `<list.txt>`: File containing one sequence name per line.
*   `-i, --invert`: Invert selection (output sequences NOT in the list).
*   `-o, --outfile <file>`: Output filename (default: stdout).
*   Output is FASTA with one sequence line per record.

---

## Transform Commands

### to-fa

Converts a 2bit file to FASTA format.

```bash
pgr 2bit to-fa [OPTIONS] <infile>
```

*   `-l, --line <int>`: Sequence line length (default: 60). Set to 0 for single line.
*   `--no-mask`: Convert sequence to all uppercase (remove soft-masking).
*   `-o, --outfile <file>`: Output filename (default: stdout).
*   Masking is preserved by default; use `--no-mask` to emit uppercase only.
