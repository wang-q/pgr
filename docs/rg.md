# pgr rg

`pgr rg` provides line-oriented operations on **`.rg` range files**
(`chr:start-end`, 1-based inclusive; species/strand prefixes are dropped).
It was migrated from `pgr runlist cover/coverage`, which in turn came from
the external `spanr` tool. Set-level runlist JSON operations live under
`pgr runlist`.

Lines starting with `#` are treated as comments and skipped by every
subcommand; lines without a valid range are skipped (or, for `sort`, written
to the end of the output).

## Subcommands

*   `cover`: Merge `.rg` range lines into a runlist JSON (per-chromosome union).
*   `coverage`: Compute depth of coverage over `.rg` ranges (sweep-line, O(n log n)).
*   `count`: Count, for each range, the overlaps with other `.rg` range files.
*   `merge`: Cluster nearly-identical `.rg` ranges and emit mappings.
*   `prop`: Proportion of each range covered by a runlist.
*   `runlist`: Filter `.rg` lines by comparing with a runlist file.
*   `sort`: Sort `.rg` lines by chromosome, start and strand.
*   `span`: Apply line-level span operations (trim/pad/shift/flank/excise).

---

## cover

Merges `chr:start-end` lines from one or more `.rg` files into a runlist JSON.
Species/strand prefixes in the lines are dropped (e.g. `S288c.I(-):190-200`
contributes to chromosome `I`).

```bash
pgr rg cover [OPTIONS] <infiles>...
pgr rg cover a.rg b.rg -o out.json
```

## coverage

Computes per-position coverage depth over `.rg` ranges with a sweep line over
sorted start/end events (O(n log n); no interval tree needed for pure depth).
Writes regions whose depth reaches `--minimum`. With `--detailed` the regions
are grouped by their exact depth instead.

```bash
pgr rg coverage [OPTIONS] <infiles>...
pgr rg coverage a.rg -m 4 -o cov.json
pgr rg coverage a.rg -m 2 -d -o cov.json
```

## count

Counts, for each range in the target `.rg` file, how many ranges in the other
`.rg` files overlap it, appending the count as an extra tab-separated field.
Lines without a valid range are skipped.

```bash
pgr rg count <target.rg> <infiles>...
pgr rg count target.rg intervals.rg
pgr rg count target.rg stdin
```

## prop

For each range in the `.rg` files, appends the proportion of the range covered
by the runlist (intersection size / range length, 4 decimals). With `--full`
the range length and the intersection size are appended as well. Lines
without a valid range are skipped.

```bash
pgr rg prop <runlist.json> <infiles>...
pgr rg prop intergenic.json a.rg
pgr rg prop intergenic.json a.rg --full
```

## sort

Sorts `.rg` lines by the parsed (chromosome, start, strand) key. Lines without
a valid range are written to the end of the output in their original order;
lines with identical keys keep their input order (stable sort).

```bash
pgr rg sort <infiles>...
pgr rg sort a.rg
pgr rg sort a.rg b.rg -o sorted.rg
```

## runlist

Keeps `.rg` lines whose range overlaps, does not overlap, or is fully contained
by the runlist, according to `--op`. `--op superset` keeps only lines whose
range is entirely inside the runlist (the runlist is a superset of the range).
Lines without a valid range are skipped.

```bash
pgr rg runlist <runlist.json> <infiles>...
pgr rg runlist intergenic.json a.rg
pgr rg runlist intergenic.json a.rg --op non-overlap
pgr rg runlist intergenic.json a.rg --op superset
```

## span

Applies an operation to each `.rg` line and writes the new range (or appends
it with `--append`): `trim`/`pad` remove or add N bases at both/5p/3p ends,
`shift`/`flank` work on the 5p or 3p end, `excise` drops ranges smaller than
N (written as an empty line).

```bash
pgr rg span <infiles>... [--op OP] [-m MODE] [-n N] [-a]
pgr rg span a.rg --op trim -n 10
pgr rg span a.rg --op flank -m 3p -n=-1 -a
```

## merge

Clusters `.rg` ranges whose reciprocal overlap reaches `--coverage` and emits
`range<TAB>merged` mapping lines for ranges in multi-member clusters. The
merged representative is the union cover `chr(+):min-max`; ranges not joined
with any other are omitted. Adapted from `rgr merge` to single-column `.rg`
input.

```bash
pgr rg merge <infiles>... [--coverage FLOAT]
pgr rg merge a.rg
pgr rg merge a.rg b.rg --coverage 0.90 -o map.tsv
```
