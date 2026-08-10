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

The bbduk `stats=` text files (`R.trim.stats.txt` / `R.filter.stats.txt`) are
not committed; regenerate them with the commands above if a stats comparison
is ever needed.

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
