# pgr CLI Examples and Notes

This document collects extended examples, end-to-end workflows, and contributor notes
that do not fit into the main README. All commands assume that `pgr` has been built
and installed as described in the README.

---

## Examples

### FA files

```bash
pgr fa size tests/fasta/ufasta.fa
pgr fa count tests/fasta/ufasta.fa.gz
pgr fa masked tests/fasta/ufasta.fa
pgr fa n50 tests/fasta/ufasta.fa -N 90 -N 50 -S

pgr fa size tests/genome/mg1655.fa.gz -o tests/genome/mg1655.size.tsv

pgr fa one tests/fasta/ufasta.fa read12
pgr fa some tests/fasta/ufasta.fa tests/fasta/list.txt
pgr fa order tests/fasta/ufasta.fa tests/fasta/list.txt

pgr fa filter tests/fasta/ufasta.fa --min-len 10 --max-len 50 --uniq
pgr fa filter tests/fasta/ufasta.fa tests/fasta/ufasta.fa.gz --min-len 1 --uniq
pgr fa filter tests/fasta/filter.fa --iupac --upper

pgr fa dedup tests/fasta/dedup.fa
pgr fa dedup tests/fasta/dedup.fa --seq --both --dups-file stdout

pgr fa mask tests/fasta/ufasta.fa --runlist tests/fasta/mask.json --hard

pgr fa replace tests/fasta/ufasta.fa --replace-tsv tests/fasta/replace.tsv
pgr fa rc tests/fasta/ufasta.fa

pgr fa filter tests/fasta/ufasta.fa --min-len 400 |
    pgr fa split name stdin -o tmp
pgr fa split about tests/fasta/ufasta.fa -c 2000 -o tmp

pgr fa six-frame tests/fasta/trans.fa
pgr fa six-frame tests/fasta/trans.fa --min-len 3 --start-met --end
```

### Block FA files

```bash
pgr maf to-fas tests/maf/example.maf

pgr axt to-fas tests/axt/RM11_1a.sizes tests/axt/example.axt --q-name RM11_1a

pgr fas filter tests/fas/example.fas --min-len 10

pgr fas name tests/fas/example.fas --count

pgr fas cover tests/fas/example.fas

pgr fas cover tests/fas/example.fas --name S288c --trim 10

pgr fas concat tests/fas/example.fas -R tests/fas/name.lst

pgr fas subset tests/fas/example.fas -R tests/fas/name.lst
pgr fas subset tests/fas/refine.fas -R tests/fas/name.lst --strict

pgr fas link tests/fas/example.fas --pair
pgr fas link tests/fas/example.fas --best

pgr fas replace tests/fas/example.fas --replace-tsv tests/fas/replace.tsv
pgr fas replace tests/fas/example.fas --replace-tsv tests/fas/replace.fail.tsv

pgr fa range tests/fas/NC_000932.fa NC_000932:1-10

pgr fas check tests/fas/A_tha.pair.fas -r tests/fas/NC_000932.fa
pgr fas check tests/fas/A_tha.pair.fas --name A_tha -r tests/fas/NC_000932.fa

pgr fas create tests/fas/I.connect.tsv -g tests/fas/genome.fa --name S288c

# Create a fasta file containing multiple genomes
cat tests/fas/genome.fa | sed 's/^>/>S288c./' > tests/fas/genomes.fa
samtools faidx tests/fas/genomes.fa S288c.I:1-100

pgr fas create tests/fas/I.connect.tsv -g tests/fas/genomes.fa

pgr fas separate tests/fas/example.fas -o . --suffix .tmp

spoa tests/fas/refine.fasta -r 1

pgr fas consensus tests/fas/example.fas
pgr fas consensus tests/fas/refine.fas
pgr fas consensus tests/fas/refine.fas --outgroup -p 2

pgr fas refine tests/fas/example.fas
pgr fas refine tests/fas/example.fas --engine none --chop 10
pgr fas refine tests/fas/refine2.fas --engine clustalw --outgroup
pgr fas refine tests/fas/example.fas --quick

pgr fas split tests/fas/example.fas --simple
pgr fas split tests/fas/example.fas -o . --chr --suffix .tmp

pgr fas slice tests/fas/slice.fas --runlist tests/fas/slice.json --name S288c

pgr fas join tests/fas/S288cvsYJM789.slice.fas --name YJM789
pgr fas join \
    tests/fas/S288cvsRM11_1a.slice.fas \
    tests/fas/S288cvsYJM789.slice.fas \
    tests/fas/S288cvsSpar.slice.fas

pgr fas stat tests/fas/example.fas --outgroup

pgr fas variation tests/fas/example.fas
pgr fas variation tests/fas/example.fas --outgroup

# snp-sites -v tests/fas_vcf/YDL184C.fas
pgr fas to-vcf tests/fas_vcf/YDL184C.fas
pgr fas to-vcf tests/fas/example.fas
pgr fas to-vcf --sizes tests/fas_vcf/S288c.chr.sizes tests/fas_vcf/YDL184C.fas

#fasops xlsx tests/fas/example.fas -o example.xlsx
#fasops xlsx tests/fas/example.fas -l 50 --outgroup -o example.outgroup.xlsx
pgr fas to-xlsx tests/fas/example.fas --indel
pgr fas to-xlsx tests/fas/example.fas --indel --outgroup
pgr fas to-xlsx tests/fas/example.fas --no-single
pgr fas to-xlsx tests/fas/example.fas --indel --no-complex
pgr fas to-xlsx tests/fas/example.fas --indel --min 0.3 --max 0.7

pgr pl p2m tests/fas/S288cvsRM11_1a.slice.fas tests/fas/S288cvsSpar.slice.fas
```

### 2bit

```bash
# pgr fa to-2bit tests/fasta/ufasta.fa -o tests/fasta/ufasta.2bit
faToTwoBit tests/genome/mg1655.fa.gz tests/genome/mg1655.2bit

pgr 2bit size tests/genome/mg1655.2bit
pgr 2bit size tests/genome/mg1655.2bit --no-ns

pgr 2bit to-fa tests/genome/mg1655.2bit -o tests/genome/mg1655.fa
pgr 2bit to-fa tests/genome/mg1655.2bit --no-mask -o tests/genome/mg1655.unmasked.fa

pgr 2bit range tests/genome/mg1655.2bit NC_000913:1-100
pgr 2bit range tests/genome/mg1655.2bit NC_000913(-):1-100

pgr 2bit masked tests/genome/mg1655.2bit
pgr 2bit masked tests/genome/mg1655.2bit --gap
```

### pbit

```bash
# Create a pbit archive from a reference and one or more samples
pgr pbit create -r ref.fa -i sample1.fa -i sample2.fa -o out.pbit

# Append samples to an existing archive
pgr pbit append out.pbit -i sample3.fa

# Show archive overview
pgr pbit stat out.pbit

# List all samples or reference contigs
pgr pbit stat out.pbit --samples
pgr pbit stat out.pbit --refs

# Extract a region from all samples
pgr pbit range out.pbit "chr1:1-1000" -o region.fa

# Extract specific contigs from all samples
pgr pbit some out.pbit contigs.txt -o selected.fa

# Export all samples as per-sample FASTA files
pgr pbit to-fa out.pbit -o outdir

# Append a new reference genome to an existing archive
pgr pbit append-ref out.pbit -r ref2.fa -o out2.pbit

# CIGAR-driven encoding with PAF (recommended for samples with SVs)
minimap2 -cx asm20 --eqx ref.fa sample.fa > sample.paf
pgr pbit create -r ref.fa -i sample.fa -p sample.paf -o out.pbit
```

### Genomes and plots

#### Distance with pgr

```bash
pgr dist mini tests/pgr/sakai.fa.gz tests/pgr/mg1655.fa.gz --hasher mod -k 21 -w 1
#NC_002695       NC_000913       0.0221  0.4580  0.5881
#NC_002127       NC_000913       0.6640  0.0000  0.0006
#NC_002128       NC_000913       0.4031  0.0001  0.0053

pgr fa rc tests/pgr/mg1655.fa.gz |
    pgr dist mini tests/pgr/sakai.fa.gz stdin --hasher mod -k 21 -w 1
#NC_002695       RC_NC_000913    0.0221  0.4580  0.5881
#NC_002127       RC_NC_000913    0.6640  0.0000  0.0006
#NC_002128       RC_NC_000913    0.4031  0.0001  0.0053

pgr fa rc tests/pgr/mg1655.fa.gz |
    pgr dist mini tests/pgr/mg1655.fa.gz stdin --hasher mod -k 21 -w 1
#NC_000913       RC_NC_000913    0.0000  1.0000  1.0000
pgr fa rc tests/pgr/mg1655.fa.gz |
    pgr dist mini tests/pgr/mg1655.fa.gz stdin --hasher rapid -k 21 -w 1
#NC_000913       RC_NC_000913    0.2289  0.0041  0.0082

pgr dist mini tests/pgr/sakai.fa.gz tests/pgr/mg1655.fa.gz --merge --hasher mod -k 21 -w 1
#tests/pgr/sakai.fa.gz   tests/pgr/mg1655.fa.gz  5302382 4543891 3064483 6781790 0.0226  0.4519  0.5779

pgr dist mini tests/pgr/sakai.fa.gz tests/pgr/mg1655.fa.gz --merge --hasher rapid -k 21 -w 1
#tests/pgr/sakai.fa.gz   tests/pgr/mg1655.fa.gz  5394043 4562542 3071076 6885509 0.0230  0.4460  0.5693

echo -e "tests/pgr/sakai.fa.gz\ntests/pgr/mg1655.fa.gz" |
    pgr dist mini stdin --merge --list-files --hasher mod -k 21 -w 1
#tests/pgr/sakai.fa.gz   tests/pgr/sakai.fa.gz   5302382 5302382 5302382 5302382 0.0000  1.0000  1.0000
#tests/pgr/sakai.fa.gz   tests/pgr/mg1655.fa.gz  5302382 4543891 3064483 6781790 0.0226  0.4519  0.5779
#tests/pgr/mg1655.fa.gz  tests/pgr/sakai.fa.gz   4543891 5302382 3064483 6781790 0.0226  0.4519  0.6744
#tests/pgr/mg1655.fa.gz  tests/pgr/mg1655.fa.gz  4543891 4543891 4543891 4543891 0.0000  1.0000  1.0000

```

### K-mer analysis (pgr kmer)

See [kmer.md](kmer.md) for the full command reference. Quick examples:

```bash
# Count k-mers into a reusable table
pgr kmer table reads.fq.gz -k 21 -o reads.pkt

# Frequency histogram (FastK .hist layout, readable by Histex/GenomeScope)
pgr kmer hist -t reads.pkt -o reads.hist
Histex -G reads.hist

# GC-content x coverage matrix (KatGC .kgc format) or its LaTeX heatmap
pgr kmer gc -t reads.pkt -o reads.kgc
pgr kmer gc -t reads.pkt --tex -o reads.tex
pgr plot heat reads.kgc -o heat.tex

# Quality-weighted histogram and error-read filtering (quorum semantics)
pgr kmer qhist reads.fq.gz -k 21 -o reads.qhist
pgr fq s-filter reads.fq.gz -k 21 -o kept.fq --discard-file bad.fq

# Coverage peak and genome-size estimate
pgr kmer gsize reads.fq.gz -k 21 --model --plot -o gs_out
pgr plot spectra reads.hist gs_out/model.txt -o spectra.tex

# Compile the figures with tectonic
tectonic heat.tex spectra.tex
```

### Alignment pipelines (lastz, UCSC, multiz)

Repeat masking (libraries, RepeatMasker, and pgr's `rept e-kmer/s-kmer/trf`) is
documented in [rept.md](rept.md).

```bash
lastz tests/pgr/pseudocat.fa tests/pgr/pseudopig.fa |
    lavToPsl stdin stdout \
    > tests/pgr/lastz.psl

pgr pl ucsc tests/pgr/pseudocat.fa tests/pgr/pseudopig.fa tests/pgr/lastz.psl

lastz --self <(gzip -dcf tests/pgr/mg1655.fa.gz)

multiz M=10 tests/multiz/S288cvsRM11_1a.maf tests/multiz/S288cvsSpar.maf 1 out1 out2
```

### Proteomes and hypervectors

```bash
pgr dist hv tests/clust/IBPA.fa
pgr dist mini tests/clust/IBPA.fa --merge

pgr dist hv tests/genome/mg1655.pro.fa.gz
pgr dist mini tests/genome/mg1655.pro.fa.gz --merge

pgr dist hv tests/genome/mg1655.pro.fa.gz tests/genome/pao1.pro.fa.gz -k 7 -w 1
pgr dist mini tests/genome/mg1655.pro.fa.gz tests/genome/pao1.pro.fa.gz -k 7 -w 1 --merge
```

### Assemblies

```bash
pgr pl prefilter tests/index/final.contigs.fa tests/clust/IBPA.fa

# SRR6323163 - APH(3')-IIIa
# 3300030246 - acrB
pgr pl prefilter tests/metagenome/SRR6323163.fa.gz "tests/metagenome/APH(3')-IIIa.fa"
pgr pl prefilter tests/metagenome/SRR6323163.fa.gz "tests/metagenome/acrB.fa"

pgr fa range tests/metagenome/SRR6323163.fa.gz "k141_4576(-):285-455|frame=2"

pgr pl prefilter tests/metagenome/3300030246.fna.gz "tests/metagenome/APH(3')-IIIa.fa" -c 1000000 -p 8
```

---

### SD (segmental duplication)

```bash
# Full SD pipeline (pgi engine, fully native):
# search -> align -> cluster -> decompose -> cover
pgr sd run tests/genome/mg1655.fa.gz -o sd_out/

# Step-by-step with intermediate files
pgr sd search tests/genome/mg1655.fa.gz --engine pgi -o hits.psl
pgr sd align tests/genome/mg1655.fa.gz hits.psl -o hits.paf
pgr sd cluster tests/genome/mg1655.fa.gz hits.paf -o clusters/
pgr sd decompose clusters/cluster_1.fa -o cluster_1.elem.bed
pgr sd cover hits.paf elems.bed -o covered.bed

# Cross-genome SD mapping
pgr sd cross tests/genome/mg1655.fa.gz tests/genome/sakai.fa.gz --engine pgi -o cross.paf
```

See `docs/sd.md` for the full command reference.

## External Dependencies (Details)

Some subcommands in `pgr` invoke external executables. In particular:

- `pgr pl ucsc` depends on the UCSC kent-tools suite, including:
  - `faToTwoBit`, `axtChain`, `chainAntiRepeat`, `chainMergeSort`, `chainPreNet`,
    `chainNet`, `netSyntenic`, `netChainSubset`, `chainStitchId`, `netSplit`,
    `netToAxt`, `axtSort`, `axtToMaf`, `netFilter`, `netClass`, `chainSplit`.
- `pgr rept trf` depends on `trf`.
- `pgr fas refine` depends on a multiple sequence alignment tool:
  - `clustalw` (default), or `muscle`, or `mafft`.

Ensure these tools are installed and available in your `PATH` before running the corresponding pipelines.
