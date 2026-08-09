# K-mer file formats (pgr kmer)

The `pgr kmer` commands produce four file formats:

* `.pkt` — canonical k-mer count table (pgr-native, single file).
* `.pkp` — per-sequence k-mer profiles (pgr-native, single file).
* `.hist` — k-mer frequency histogram (FastK-compatible binary).
* `.kgc` — GC-content × coverage matrix (KatGC-compatible TSV).

## `.pkt` (k-mer table)

Native format of `pgr kmer table`, also used as the `rept e-kmer`
repeat-library cache (magic `PKTT`). Layout:

```text
header (bincode, 24 bytes):
  magic      [u8; 4]   = "PKTT"
  version    u32       = 1
  k          u32       (k-mer length)
  n_entries  u64
  key_bytes  u32       = ceil(2k / 8)
entry, repeated n_entries times:
  packed key key_bytes bytes (big-endian 2-bit k-mer, low bits zero-padded)
  count      u32 (little-endian)
```

Keys are canonical (the lexicographically smaller of the forward and
reverse-complement 2-bit encodings), ascending and duplicate-free. Written
atomically (temp file + rename). Read with `libs/kmer/count.rs`
(`save`/`load`/`k_of`).

## `.pkp` (per-sequence profile)

Native format of `pgr kmer profile` (magic `PKPP`). Layout:

```text
header:
  magic      [u8; 4]   = "PKPP"
  version    u32       = 1
  k          u32
  n_seqs     u64
per sequence:
  length     u64       (profile length = seq len - k + 1)
  values     u16[length] (little-endian, raw counts)
```

`values[i]` is the count of the canonical k-mer at position `i` (0 for N
windows); self profiles use dataset-wide counts, relative profiles use the
supplied table's counts. Read with `libs/kmer/profile.rs`
(`save_profiles`/`load_profiles`).

## `.hist` (frequency histogram)

Byte-compatible with the FastK `.hist` binary layout (produced by
`pgr kmer hist`; read by Histex / KatGC / GenomeScope tooling). Fixed
28-byte header plus 32767 `int64` bins (262,164 bytes total):

```text
kmer      int32   (k-mer length)
low       int32   = 1
high      int32   = 32767
ilowcnt   int64   (= bin 1 count)
max_inst  int64   (instance-mode top-bin boundary)
hist      int64[32767]   hist[i-1] = distinct k-mers with count i
```

All integers little-endian. Counts above 32767 fold into the top bin and
their instances accumulate in `max_inst`, matching FastK. Implemented in
`libs/kmer/hist.rs`; verified identical to a real FastK histogram.

## `.kgc` (GC × coverage matrix)

TSV output of `pgr kmer gc`, matching MerquryFK `KatGC`'s temporary
`.kgc` matrix:

```text
GCP\tKF\tCount
```

* `GCP` — GC-count bin center (`i.5`, `i` = GC count of the k-mer).
* `KF` — coverage bin center (`a.5`, `a` = count bin).
* `Count` — 2×2 neighbor average of the matrix, clamped to the peak value.

Rows span GC counts `0..k-1`, columns `0..xmax-1` (peak × `--xrel`, default
2.1). Implemented in `libs/kmer/gc.rs`; verified line-identical to a
locally compiled KatGC.
