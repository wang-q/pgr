# pgr net

`pgr net` provides tools for manipulating **Net** alignment files. Net files (UCSC format) represent the hierarchical structure of alignments.

## Subcommands

*   `class`: Show statistics of net classes (e.g., top, syn, nonSyn).
*   `filter`: Filter net files based on score, size, and synteny criteria.
*   `split`: Split a net file into individual files per chromosome.
*   `subset`: Create a chain file containing only the chains referenced in the net.
*   `syntenic`: Add synteny information (class labels) to a net file.
*   `to-axt`: Convert net and chain files to AXT format.

---

## class

Displays statistics of the net file, including the count and total bases for each class (e.g., `top`, `syn`, `inv`, `nonSyn`, `gap`).

```bash
pgr net class [input]
```

### Options

*   `[input]`: Input net file (or `-` for stdin).

### Examples

1.  **Show net statistics**:
    ```bash
    pgr net class in.net
    ```

---

## filter

Filters a net file based on various criteria such as score, gap size, alignment size, and sequence names. It can also filter based on synteny status.

```bash
pgr net filter [OPTIONS] <input>
```

### Options

*   `--min-score <float>`: Restrict to fills scoring at least this value.
*   `--max-score <float>`: Restrict to fills scoring at most this value.
*   `--min-gap <int>`: Restrict to gaps with size >= this value.
*   `--min-ali <int>`: Restrict to fills with at least this many aligned bases.
*   `--max-ali <int>`: Restrict to fills with at most this many aligned bases.
*   `--min-size-t <int>`: Restrict to fills with target size >= this value.
*   `--min-size-q <int>`: Restrict to fills with query size >= this value.
*   `--t <names>`: Restrict target sequence to those named (comma-separated).
*   `--not-t <names>`: Restrict target sequence to those NOT named (comma-separated).
*   `--q <names>`: Restrict query sequence to those named (comma-separated).
*   `--not-q <names>`: Restrict query sequence to those NOT named (comma-separated).
*   `--type <type>`: Restrict to given type (e.g., `top`, `syn`). Can be repeated.
*   `--syn`: Filter based on synteny criteria (tuned for human/mouse).
*   `--nonsyn`: Inverse filtering based on synteny.
*   `--fill-only`: Only pass fills, not gaps.
*   `--gap-only`: Only pass gaps, not fills.

### Examples

1.  **Filter for high-scoring syntenic blocks**:
    ```bash
    pgr net filter in.net --syn --min-score 10000 > out.net
    ```

---

## split

Splits a single net file containing multiple chromosomes into individual files (one per chromosome) in the specified directory.

```bash
pgr net split <input> <output_dir>
```

### Examples

1.  **Split net file by chromosome**:
    ```bash
    pgr net split all.net nets/
    ```

---

## subset

Creates a new chain file that contains only the chains (or parts of chains) referenced by the net file. This is useful for reducing the size of chain files to only the "best" alignments selected by the net.

```bash
pgr net subset [OPTIONS] <net_in> <chain_in> <chain_out>
```

### Options

*   `--whole-chains`: Write entire chains referenced by the net (don't split/subset).
*   `--split-on-insert`: Split chain when an insertion of another chain occurs (nested structure).
*   `--type <string>`: Restrict output to chains associated with a specific net type.

### Examples

1.  **Create a subset chain file**:
    ```bash
    pgr net subset in.net in.chain out.chain
    ```

---

## syntenic

Adds synteny information to a net file. It classifies fills as `syn` (syntenic), `inv` (inverted), or `nonSyn` (non-syntenic) based on their relationship to the parent structure and query coordinates.

```bash
pgr net syntenic [OPTIONS] <in_net> <out_net>
```

### Options

*   `--min-score <float>`: Minimum score to output (default: 0.0).

### Examples

1.  **Add synteny labels**:
    ```bash
    pgr net syntenic raw.net syn.net
    ```

---

## to-axt

Converts a net file (and its corresponding chain file) into AXT format. This reconstructs the alignment sequence from the 2bit files.

```bash
pgr net to-axt <in_net> <in_chain> <target.2bit> <query.2bit> <out_axt>
```

### Arguments

*   `<in_net>`: Input net file.
*   `<in_chain>`: Input chain file.
*   `<target.2bit>`: Target sequence 2bit file.
*   `<query.2bit>`: Query sequence 2bit file.
*   `<out_axt>`: Output AXT file.

### Examples

1.  **Convert net to AXT**:
    ```bash
    pgr net to-axt in.net in.chain target.2bit query.2bit out.axt
    ```
