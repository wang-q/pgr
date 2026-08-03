# pgr runlist

`pgr runlist` provides interval-set operations on **runlists** — either `.rg`
lines (`chr:start-end`, 1-based inclusive) or runlist JSON
(`{"chr": "start-end,..."}`), the format consumed by `pgr fa mask`. The
command family was migrated from the external `spanr` tool (intspan project)
and produces identical output.

## Subcommands

*   `cover`: Merge `.rg` range lines into a runlist JSON (per-chromosome union).
*   `coverage`: Compute depth of coverage over `.rg` ranges (sweep-line, O(n log n)).
*   `span`: Apply span operations (cover/holes/trim/pad/excise/fill) to a runlist JSON.
*   `compare`: Set operations (intersect/union/diff/xor) between runlist JSON files.
*   `merge`: Merge runlist JSON files into a multi runlist keyed by file stem.

---

## cover

Merges `chr:start-end` lines from one or more `.rg` files into a runlist JSON.
Species/strand prefixes in the lines are dropped (e.g. `S288c.I(-):190-200`
contributes to chromosome `I`).

```bash
pgr runlist cover [OPTIONS] <infiles>...
pgr runlist cover a.rg b.rg -o out.json
```

## coverage

Computes per-position coverage depth over `.rg` ranges with a sweep line over
sorted start/end events (O(n log n); no interval tree needed for pure depth).
Writes regions whose depth reaches `--minimum`. With `--detailed` the regions
are grouped by their exact depth instead.

```bash
pgr runlist coverage [OPTIONS] <infiles>...
pgr runlist coverage a.rg -m 4 -o cov.json
pgr runlist coverage a.rg -m 2 -d -o cov.json
```

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
