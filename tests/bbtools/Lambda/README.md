# Lambda golden data for the BBTools migration

Source: `anchr/tests/Lambda/` (reads SRR5042715, lambda phage Illumina PE).

`R1.fq.gz` / `R2.fq.gz` are the raw inputs; `illumina_adapters.fa` is the
reference file used by the trim/filter steps of the pipeline.

`golden/` holds the byte-level reference outputs of every BBTools step, produced
with the locally installed **BBTools 39.38** (cbp package, entry `jgi.BBDuk`)
using the anchr `trim.era.sh` defaults:

* `--trimk 23 --matchk 27 --cutk 31 --trimq 15 --qual 25 --len 60`
* `--filter adapter`, no dedupe / tile / cutoff / sampling

Determinism notes (both required for reproducible output):

* bbduk runs use `ordered=t`; the BBTools default is unordered and
  nondeterministic across runs.
* clumpify and reformat use fixed seeds (`seed=1`, `sampleseed=1`).

Generation commands (run in one work directory, `in/` holding the inputs):

```bash
clumpify.sh  -Xmx2g in=in/R1.fq.gz in2=in/R2.fq.gz out=clumpify.fq.gz \
    threads=8 seed=1
bbduk.sh -Xmx2g in=clumpify.fq.gz out=trim.fq.gz \
    ref=in/illumina_adapters.fa \
    maxns=0 ktrim=r k=23 mink=11 hdist=1 tbo tpe \
    minlen=60 qtrim=r trimq=15 ftm=5 \
    stats=R.trim.stats.txt overwrite tossbrokenreads=t threads=8 ordered=t
bbduk.sh -Xmx2g in=trim.fq.gz out=filter.fq.gz \
    ref=in/illumina_adapters.fa k=27 cardinality \
    stats=R.filter.stats.txt overwrite tossbrokenreads=t threads=8 ordered=t
kmercountexact.sh -Xmx2g in=filter.fq.gz khist=R.khist.txt \
    peaks=R.peaks.txt k=31 threads=8
repair.sh -Xmx2g in=filter.fq.gz out=out/R1.fq.gz out2=out/R2.fq.gz \
    outs=out/Rs.fq.gz repair threads=8
reformat.sh -Xmx2g in=filter.fq.gz out=sample.fq.gz \
    samplebasestarget=1000000 sampleseed=1 threads=8
```

All FASTQ goldens are stored gzipped (repair outputs `R1.fq.gz` / `R2.fq.gz` /
`Rs.fq.gz` were re-gzipped from the repair outputs). Comparisons in pgr tests
decompress both sides, since gzip bytes (mtime) are not stable.

## merge (bbmerge)

`anchr fq merge` goldens (`merge.*.fq.gz`, `merge.*.txt`) were produced with the
locally installed **BBTools 40.01** (`jgi.BBMerge`) using `ordered=t threads=1`
on a 2000-pair subset (`R1.2k.fq.gz` / `R2.2k.fq.gz`, the first 2000 pairs of
`R1.fq.gz` / `R2.fq.gz`). bbmerge processes each pair independently, so the
2k-subset goldens are byte-identical to the corresponding records of a
40000-pair run; the subset keeps the test data ~10x smaller
(the merge flow needs the bundled `bbmerge.bbnet`, committed here as
`bbmerge.bbnet`):

```bash
# net path (vstrict ecco, anchr merge phase 1; default makevector=true)
bbmerge.sh in=R1.2k.fq in2=R2.2k.fq out=merge.ecco.fq.gz \
    ihist=merge.ihist1.txt threads=1 ecco mix vstrict ordered=t overwrite

# net path (strict join)
bbmerge.sh in=R1.2k.fq in2=R2.2k.fq out=merge.merged.fq.gz outu=merge.unmerged.fq.gz \
    ihist=merge.ihist2.txt threads=1 strict ordered=t overwrite

# classic path (makevector=f, no net)
bbmerge.sh in=R1.2k.fq in2=R2.2k.fq out=merge.novector.merged.fq.gz \
    outu=merge.novector.unmerged.fq.gz ihist=merge.novector.ihist.txt \
    threads=1 strict makevector=f ordered=t overwrite
bbmerge.sh in=R1.2k.fq in2=R2.2k.fq out=merge.novector.ecco.fq.gz \
    threads=1 ecco mix vstrict makevector=f ordered=t overwrite
```

The `bbmerge.sh` wrapper adds `-Xmx` and classpath flags; the equivalent direct
invocation is `java -Xmx3g -cp BBTools-40.01/current jgi.BBMerge ...`. The
ecco output is multi-member gzip when written by BBTools, so the committed
`merge.ecco.fq.gz` was re-gzipped to a single member (as with the repair
outputs); pgr tests decompress with `MultiGzDecoder` to be safe.

Note: BBTools' paired-file reader (`PairStreamer`) desynchronizes batch sizes
on gzipped `in=`/`in2=` inputs, so the golden commands above use plain-text
`R1.2k.fq` / `R2.2k.fq` (uncompressed copies of the subset) and the outputs
are gzipped afterwards; pgr reads the gzipped inputs directly.

The bbduk `stats=` text files (`R.trim.stats.txt` / `R.filter.stats.txt`) are
not committed; regenerate them with the commands above if a stats comparison
is ever needed. `anchr fq clean --stats` reproduces the 3-column format
byte for byte (values checked in `cli_fq_trim_adapter.rs`); the `#File` line
carries the input path, so it is path-dependent by design.

## bbnorm cutoff (optional pipeline step)

The trim pipeline only runs `bbnorm.sh` when a cutoff is configured (off by
default). It is covered by the small synthetic `kmer norm` test only; the
full-data commands are kept here for reference:

```bash
bbnorm.sh -Xmx2g in=clumpify.fq.gz out=norm_min3.fq.gz \
    passes=1 bits=16 min=3 target=9999999 threads=8 ordered=t
```

Neither the BBTools output nor a pgr regression golden is committed.
`pgr kmer norm` uses an exact canonical count table instead of bbnorm's
`bits=16` approximate hash counts, so reads at the depth-3 boundary can
differ from `bbnorm.sh` output.

## merge pipeline: tadpole ecc and extend (golden committed)

`anchr fq ecc` and `anchr fq extend` were byte-compared against BBTools 40.01 on
the full Lambda data (40000 pairs), then committed as a 2000-pair subset
golden (`ecco_sub.fq.gz` input, `ecct_sub.fq.gz` / `ext_sub.fq.gz` outputs):

```bash
# Error correction + tossing (anchr merge phase 3)
java -Xmx4g -cp BBTools-40.01/current assemble.Tadpole \
    in=ecco_sub.fq out=ecct_sub.fq threads=1 \
    ecc tossjunk tossdepth=2 tossuncorrectable overwrite

# Read extension (anchr "Read extension" step, k=62 -> Tadpole2 path)
java -Xmx4g -cp BBTools-40.01/current assemble.Tadpole \
    in=ecco_sub.fq out=ext_sub.fq threads=1 \
    mode=extend el=20 er=20 k=62 overwrite
```

Both commands match the full 40000-pair run byte for byte (all reads,
sequences, qualities, and discard decisions). Notable semantics reproduced:
N bases reset the k-mer window and the minprob product everywhere (the table
never contains N-spanning windows); absent k-mers read as -1 in the count
arrays; error correction compares the base code to the *count*
(`if(num==rightMax)`) at the reassembly step; read extension never uses the
left-counts junction check (`leftCounts` is null in BBTools), and the
junction-base append condition flips between Tadpole1 (`kmer>rkmer`) and
Tadpole2 (`kmer<rkmer`).

## merge phase 4: bbmerge-auto extend2/rem (golden committed)

`anchr fq merge --extend2 N --rem` reproduces anchr merge phase 4
(`bbmerge-auto.sh ... strict k=81 extend2=80 rem`): pairs that fail the
classic overlap are extended along the k=81 k-mer graph and re-checked, with
`requireExtensionMatch` requiring the extended overlap to agree with the
unextended one. Golden (`merge4.*`) was produced with BBTools 40.01 over the
2000-pair extended reads (`ext_sub.fq.gz`):

```bash
java -Xmx4g -cp BBTools-40.01/current jgi.BBMerge \
    in=ext_sub.fq.gz out=merge4.merged.fq.gz outu=merge4.unmerged.fq.gz \
    ihist=merge4.ihist.txt threads=1 strict k=81 extend2=80 rem overwrite
```

The output is byte-identical to `anchr fq merge ... --strict --no-make-vector
--extend2 80 --rem` (merged, unmerged, and ihist). Key semantics: BBMerge
snapshots the pre-extension reads and restores them for unmerged output;
`lengthSum` for the rem acceptance rule is the *unextended* length; the
extend2 extension uses `includeJunctionBase=false` and never checks left
junctions (`extendThroughLeftJunctions` defaults true -> leftCounts null);
`extendIterations` defaults to 1; and every pair is extended when `rem` is
set, not just ambiguous ones.

## clumpify dedupe (verified once, no golden)

`anchr fq clump --dedupe --dupesubs 0` was byte-compared against BBTools 39.38
`clumpify.sh ... threads=1 dedupe=t dupesubs=0` on the Lambda data (40000 ->
39984 reads; R1 and R2 both exact, N wildcard; higher-quality copy kept) and
matched byte for byte. No golden is committed; a small synthetic test covers
the semantics. Note: with `threads>1`, BBTools dedupe output order is
nondeterministic (clump/thread race; removed set identical), so byte
comparison uses `threads=1`.
