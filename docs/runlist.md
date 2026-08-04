# pgr runlist

`pgr runlist` provides interval-set operations on **runlist JSON**
(`{"chr": "start-end,..."}`), the format consumed by `pgr fa mask`. The
command family was migrated from the external `spanr` tool (intspan project)
and produces identical output. Line-oriented `.rg` input is handled by
`pgr rg` instead.

## Subcommands

*   `combine`: Combine multiple sets of a multi runlist JSON into one.
*   `compare`: Set operations (intersect/union/diff/xor) between runlist JSON files.
*   `convert`: Convert runlist JSON files to `.rg` range lines.
*   `genome`: Convert a chromosome sizes file to a full-genome runlist JSON.
*   `merge`: Merge runlist JSON files into a multi runlist keyed by file stem.
*   `some`: Extract selected top-level keys from a runlist JSON.
*   `span`: Apply span operations (cover/holes/trim/pad/excise/fill) to a runlist JSON.
*   `split`: Split a multi runlist JSON into per-key files.
*   `stat`: Per-chromosome coverage stats of a runlist against chromosome sizes.
*   `statop`: Cross-set coverage stats (one runlist compared to another).

## span

Applies an operation to every chromosome of a runlist JSON (single or multi):

*   `cover`: a single span from min to max
*   `holes`: all the holes in the runlist
*   `trim`: remove N integers from each end of each span
*   `pad`: add N integers to each end of each span
*   `excise`: remove all spans smaller than N
*   `fill`: fill in all holes smaller than or equal to N

```bash
pgr runlist span [OPTIONS] <infile>
pgr runlist span in.json --op fill -n 10 -o out.json
pgr runlist span in.json --op excise -n 100 -o out.json
```

`<infile>` may be `stdin`.

## compare

Applies a set operation between the first file (which may hold multiple
runlist sets) and each of the other files. Missing chromosomes are treated as
empty.

```bash
pgr runlist compare [OPTIONS] <infile> <infiles>...
pgr runlist compare a.json b.json c.json --op intersect -o out.json
```

## merge

Reads several runlist JSON files and writes a multi runlist keyed by file
stem. Without `--all` only the first dot-separated segment of the stem is
used as the key.

```bash
pgr runlist merge [OPTIONS] <infiles>...
pgr runlist merge a.json b.json -o out.json
```

## combine

Combines all sets of a multi runlist JSON into one, applying the operation
between the first set and each subsequent set.

```bash
pgr runlist combine [OPTIONS] <infile>
pgr runlist combine in.json -o out.json
pgr runlist combine in.json --op intersect -o out.json
```

## convert

Writes `chr:start-end` lines (one per span) for every chromosome of each
input runlist. With `--longest` only the longest span per chromosome is kept.

```bash
pgr runlist convert [OPTIONS] <infiles>...
pgr runlist convert in.json -o out.rg
pgr runlist convert in.json --longest -o out.rg
```

## genome

Builds a runlist JSON where every chromosome spans its full length (1..size).

```bash
pgr runlist genome <chr.sizes> -o genome.json
```

## some

Keeps only the top-level keys listed in the names file (one per line).

```bash
pgr runlist some <infile> <list> -o out.json
```

## split

Splits a multi runlist JSON into one JSON per top-level key, written to
`<outdir>/<key><suffix>` or printed line by line with `-o stdout`.

```bash
pgr runlist split <infile> -o out_dir
pgr runlist split <infile> --suffix .json -o stdout
```

## stat

Prints per-chromosome coverage as CSV (`key,chr,chrLength,size,coverage`
plus an `all` row). `--all` keeps only the whole-genome stats.

```bash
pgr runlist stat [OPTIONS] <chr.sizes> <infile>
pgr runlist stat chr.sizes in.json -o stat.csv
```

## statop

Prints CSV stats comparing `infile1` (possibly multi) against `infile2`
(single): `key,chr,chrLength,size,<base>Length,<base>Size,c1,c2,ratio`,
where `<base>` is the stem of `infile2` (or `--base`). `--all` keeps only the
whole-genome stats.

```bash
pgr runlist statop [OPTIONS] <chr.sizes> <infile1> <infile2>
pgr runlist statop chr.sizes a.json b.json -o statop.csv
```
