# pgr rg

`pgr rg` provides line-oriented operations on **`.rg` range files**
(`chr:start-end`, 1-based inclusive; species/strand prefixes are dropped).
It was migrated from `pgr runlist cover/coverage`, which in turn came from
the external `spanr` tool. Set-level runlist JSON operations live under
`pgr runlist`.

## Subcommands

*   `cover`: Merge `.rg` range lines into a runlist JSON (per-chromosome union).
*   `coverage`: Compute depth of coverage over `.rg` ranges (sweep-line, O(n log n)).
*   `count`: Count, for each range, the overlaps with other `.rg` range files.

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
