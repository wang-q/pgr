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

pgr fa filter tests/fasta/ufasta.fa -a 10 -z 50 --uniq
pgr fa filter tests/fasta/ufasta.fa tests/fasta/ufasta.fa.gz -a 1 --uniq
pgr fa filter tests/fasta/filter.fa --iupac --upper

pgr fa dedup tests/fasta/dedup.fa
pgr fa dedup tests/fasta/dedup.fa --seq --both --file stdout

pgr fa mask tests/fasta/ufasta.fa tests/fasta/mask.json --hard

pgr fa replace tests/fasta/ufasta.fa tests/fasta/replace.tsv
pgr fa rc tests/fasta/ufasta.fa

pgr fa filter tests/fasta/ufasta.fa -a 400 |
    pgr fa split name stdin -o tmp
pgr fa split about tests/fasta/ufasta.fa -c 2000 -o tmp

pgr fa six-frame tests/fasta/trans.fa
pgr fa six-frame tests/fasta/trans.fa --len 3 --start --end
```

### Block FA files

```bash
pgr maf to-fas tests/maf/example.maf

pgr axt to-fas tests/fas/RM11_1a.chr.sizes tests/fas/example.axt --qname RM11_1a

pgr fas filter tests/fas/example.fas --ge 10

pgr fas name tests/fas/example.fas --count

pgr fas cover tests/fas/example.fas

pgr fas cover tests/fas/example.fas --name S288c --trim 10

pgr fas concat tests/fas/example.fas -r tests/fas/name.lst

pgr fas subset tests/fas/example.fas -r tests/fas/name.lst
pgr fas subset tests/fas/refine.fas -r tests/fas/name.lst --strict

pgr fas link tests/fas/example.fas --pair
pgr fas link tests/fas/example.fas --best

pgr fas replace tests/fas/example.fas -r tests/fas/replace.tsv
pgr fas replace tests/fas/example.fas -r tests/fas/replace.fail.tsv

pgr fa range tests/fas/NC_000932.fa NC_000932:1-10

pgr fas check tests/fas/A_tha.pair.fas -r tests/fas/NC_000932.fa
pgr fas check tests/fas/A_tha.pair.fas --name A_tha -r tests/fas/NC_000932.fa

pgr fas create tests/fas/I.connect.tsv -r tests/fas/genome.fa --name S288c

# Create a fasta file containing multiple genomes
cat tests/fas/genome.fa | sed 's/^>/>S288c./' > tests/fas/genomes.fa
samtools faidx tests/fas/genomes.fa S288c.I:1-100

cargo run --bin pgr -- fas create tests/fas/I.name.tsv -r tests/fas/genomes.fa --multi

pgr fas separate tests/fas/example.fas -o . --suffix .tmp

spoa tests/fas/refine.fasta -r 1

pgr fas consensus tests/fas/example.fas
pgr fas consensus tests/fas/refine.fas
pgr fas consensus tests/fas/refine.fas --outgroup -p 2

pgr fas refine tests/fas/example.fas
pgr fas refine tests/fas/example.fas --msa none --chop 10
pgr fas refine tests/fas/refine2.fas --msa clustalw --outgroup
pgr fas refine tests/fas/example.fas --quick

pgr fas split tests/fas/example.fas --simple
pgr fas split tests/fas/example.fas -o . --chr --suffix .tmp

pgr fas slice tests/fas/slice.fas -r tests/fas/slice.json --name S288c

cargo run --bin pgr -- fas join tests/fas/S288cvsYJM789.slice.fas --name YJM789
cargo run --bin pgr -- fas join \
    tests/fas/S288cvsRM11_1a.slice.fas \
    tests/fas/S288cvsYJM789.slice.fas \
    tests/fas/S288cvsSpar.slice.fas

cargo run --bin pgr -- fas stat tests/fas/example.fas --outgroup

cargo run --bin pgr -- fas variation tests/fas/example.fas
cargo run --bin pgr -- fas variation tests/fas/example.fas --outgroup

# snp-sites -v tests/fas/YDL184C.fas
cargo run --bin pgr -- fas to-vcf tests/fas/YDL184C.fas
cargo run --bin pgr -- fas to-vcf tests/fas/example.fas
cargo run --bin pgr -- fas to-vcf --sizes tests/fas/S288c.chr.sizes tests/fas/YDL184C.fas

#fasops xlsx tests/fas/example.fas -o example.xlsx
#fasops xlsx tests/fas/example.fas -l 50 --outgroup -o example.outgroup.xlsx
pgr fas to-xlsx tests/fas/example.fas --indel
pgr fas to-xlsx tests/fas/example.fas --indel --outgroup
pgr fas to-xlsx tests/fas/example.fas --nosingle
pgr fas to-xlsx tests/fas/example.fas --indel --nocomplex
pgr fas to-xlsx tests/fas/example.fas --indel --min 0.3 --max 0.7

cargo run --bin pgr -- pl p2m tests/fas/S288cvsRM11_1a.slice.fas tests/fas/S288cvsSpar.slice.fas
```

### 2bit

```bash
# pgr fa to-2bit tests/fasta/ufasta.fa -o tests/fasta/ufasta.2bit
faToTwoBit tests/genome/mg1655.fa.gz tests/genome/mg1655.2bit

pgr 2bit size tests/genome/mg1655.2bit
pgr 2bit size tests/genome/mg1655.2bit --no-ns
pgr 2bit size tests/genome/mg1655.2bit tests/genome/sakai.2bit

pgr 2bit to-fa tests/genome/mg1655.2bit -o tests/genome/mg1655.fa
pgr 2bit to-fa tests/genome/mg1655.2bit --no-mask -o tests/genome/mg1655.unmasked.fa

pgr 2bit range tests/genome/mg1655.2bit NC_000913:1-100
pgr 2bit range tests/genome/mg1655.2bit NC_000913(-):1-100
# pgr 2bit range tests/genome/mg1655.2bit --rgfile tests/genome/ranges.txt

pgr 2bit masked tests/genome/mg1655.2bit
pgr 2bit masked tests/genome/mg1655.2bit --gap
```

### Genomes and plots

#### Download genomes

```bash
curl -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz |
    gzip -dc |
    hnsm filter stdin -s |
    hnsm gz stdin -o tests/pgr/mg1655.fa

curl -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz |
    gzip -dc |
    hnsm filter stdin -s |
    hnsm gz stdin -o tests/pgr/sakai.fa
```

#### Distance with hnsm

```bash
hnsm distance tests/pgr/sakai.fa.gz tests/pgr/mg1655.fa.gz --hasher mod -k 21 -w 1
#NC_002695       NC_000913       0.0221  0.4580  0.5881
#NC_002127       NC_000913       0.6640  0.0000  0.0006
#NC_002128       NC_000913       0.4031  0.0001  0.0053

hnsm rc tests/pgr/mg1655.fa.gz |
    hnsm distance tests/pgr/sakai.fa.gz stdin --hasher mod -k 21 -w 1
#NC_002695       RC_NC_000913    0.0221  0.4580  0.5881
#NC_002127       RC_NC_000913    0.6640  0.0000  0.0006
#NC_002128       RC_NC_000913    0.4031  0.0001  0.0053

hnsm rc tests/pgr/mg1655.fa.gz |
    hnsm distance tests/pgr/mg1655.fa.gz stdin --hasher mod -k 21 -w 1
#NC_000913       RC_NC_000913    0.0000  1.0000  1.0000
hnsm rc tests/pgr/mg1655.fa.gz |
    hnsm distance tests/pgr/mg1655.fa.gz stdin --hasher rapid -k 21 -w 1
#NC_000913       RC_NC_000913    0.2289  0.0041  0.0082

hnsm distance tests/pgr/sakai.fa.gz tests/pgr/mg1655.fa.gz --merge --hasher mod -k 21 -w 1
#tests/pgr/sakai.fa.gz   tests/pgr/mg1655.fa.gz  5302382 4543891 3064483 6781790 0.0226  0.4519  0.5779

hnsm distance tests/pgr/sakai.fa.gz tests/pgr/mg1655.fa.gz --merge --hasher rapid -k 21 -w 1
#tests/pgr/sakai.fa.gz   tests/pgr/mg1655.fa.gz  5394043 4562542 3071076 6885509 0.0230  0.4460  0.5693

echo -e "tests/pgr/sakai.fa.gz\ntests/pgr/mg1655.fa.gz" |
    hnsm distance stdin --merge --list --hasher mod -k 21 -w 1
#tests/pgr/sakai.fa.gz   tests/pgr/sakai.fa.gz   5302382 5302382 5302382 5302382 0.0000  1.0000  1.0000
#tests/pgr/sakai.fa.gz   tests/pgr/mg1655.fa.gz  5302382 4543891 3064483 6781790 0.0226  0.4519  0.5779
#tests/pgr/mg1655.fa.gz  tests/pgr/sakai.fa.gz   4543891 5302382 3064483 6781790 0.0226  0.4519  0.6744
#tests/pgr/mg1655.fa.gz  tests/pgr/mg1655.fa.gz  4543891 4543891 4543891 4543891 0.0000  1.0000  1.0000

```

#### Plotting alignments

```bash
FastGA -v -pafx tests/pgr/sakai.fa.gz tests/pgr/mg1655.fa.gz > tmp.paf
FastGA -v -psl tests/pgr/sakai.fa.gz tests/pgr/mg1655.fa.gz > tmp.psl

pgr chain -t="" tests/pgr/mg1655.fa.gz tests/pgr/sakai.fa.gz tmp.psl > tmp.chain.maf
pgr chain --syn -t="" tests/pgr/mg1655.fa.gz tests/pgr/sakai.fa.gz tmp.psl > tmp.syn.maf

lastz <(gzip -dcf tests/pgr/mg1655.fa.gz) <(gzip -dcf tests/pgr/sakai.fa.gz) |
    lavToPsl stdin stdout \
    > tmp.lastz.psl
pgr chain --syn -t="" tests/pgr/mg1655.fa.gz tests/pgr/sakai.fa.gz tmp.lastz.psl > tmp.lastz.maf

wgatools dotplot -f paf tmp.paf > tmp.html
wgatools dotplot tmp.chain.maf > tmp.chain.html
wgatools dotplot tmp.syn.maf > tmp.syn.html
wgatools dotplot tmp.lastz.maf > tmp.lastz.html
```


| ![paf.png](images/paf.png) | ![chain.png](images/chain.png) |
|:--------------------------:|:------------------------------:|
|            paf             |             chain              |

| ![syn.png](images/syn.png) | ![lastz.png](images/lastz.png) |
|:--------------------------:|:------------------------------:|
|            syn             |             lastz              |


### Repeats and RepeatMasker

```bash
# TnCentral
curl -LO https://tncentral.ncc.unesp.br/api/download_blast/nc/tn_in_is

unzip -j tn_in_is 'tncentral_integrall_isfinder.fa'
gzip -9 -c 'tncentral_integrall_isfinder.fa' > tncentral.fa.gz

hnsm size tests/pgr/tncentral.fa.gz
hnsm distance tests/pgr/tncentral.fa.gz -k 17 -w 5 -p 8 |
    rgr filter stdin --ge 5:0.9

# RepBase for RepeatMasker
curl -LO https://github.com/wang-q/ubuntu/releases/download/20190906/repeatmaskerlibraries-20140131.tar.gz

tar xvfz repeatmaskerlibraries-20140131.tar.gz Libraries/RepeatMaskerLib.embl

# https://sourceforge.net/projects/readseq/
java -jar ~/bin/readseq.jar -f fa Libraries/RepeatMaskerLib.embl
mv Libraries/RepeatMaskerLib.embl.fasta repbase.fa
gzip -9 -k repbase.fa
```

```bash
singularity run ~/bin/repeatmasker_master.sif /app/RepeatMasker/RepeatMasker \
    ./genome.fa -xsmall -species "bacteria"

singularity run ~/bin/repeatmasker_master.sif /app/RepeatMasker/util/rmOutToGFF3.pl \
    ./genome.fa.out > mg1655.rm.gff

spanr gff tests/pgr/mg1655.rm.gff -o tests/pgr/mg1655.rm.json
```

```bash
pgr ir tests/pgr/tncentral.fa.gz tests/pgr/mg1655.fa.gz \
    > tests/pgr/mg1655.ir.json

spanr stat tests/pgr/mg1655.chr.sizes tests/pgr/mg1655.ir.json

pgr rept tests/pgr/mg1655.fa.gz \
    > tests/pgr/mg1655.rept.json

pgr trf tests/pgr/mg1655.fa.gz \
    > tests/pgr/mg1655.trf.json

spanr stat tests/pgr/mg1655.chr.sizes tests/pgr/mg1655.rm.json
spanr statop tests/pgr/mg1655.chr.sizes tests/pgr/mg1655.ir.json tests/pgr/mg1655.rm.json

lastz tests/pgr/pseudocat.fa tests/pgr/pseudopig.fa |
    lavToPsl stdin stdout \
    > tests/pgr/lastz.psl

pgr chain tests/pgr/pseudocat.fa tests/pgr/pseudopig.fa tests/pgr/lastz.psl

lastz --self <(gzip -dcf tests/pgr/mg1655.fa.gz)

multiz M=10 tests/multiz/S288cvsRM11_1a.maf tests/multiz/S288cvsSpar.maf 1 out1 out2
```

### Proteomes and hypervectors

```bash
pgr dist hv tests/clust/IBPA.fa
hnsm dist seq tests/clust/IBPA.fa --merge

hnsm dist hv tests/genome/mg1655.pro.fa.gz
hnsm dist seq tests/genome/mg1655.pro.fa.gz --merge

hnsm dist hv tests/genome/mg1655.pro.fa.gz tests/genome/pao1.pro.fa.gz -k 7 -w 1
hnsm dist seq tests/genome/mg1655.pro.fa.gz tests/genome/pao1.pro.fa.gz -k 7 -w 1 --merge
```

### Assemblies

```bash
cargo run --bin pgr pl prefilter tests/index/final.contigs.fa tests/clust/IBPA.fa

# SRR6323163 - APH(3')-IIIa
# 3300030246 - acrB
pgr pl prefilter tests/metagenome/SRR6323163.fa.gz "tests/metagenome/APH(3')-IIIa.fa"
pgr pl prefilter tests/metagenome/SRR6323163.fa.gz "tests/metagenome/acrB.fa"

pgr fa range tests/metagenome/SRR6323163.fa.gz "k141_4576(-):285-455|frame=2"

pgr pl prefilter tests/metagenome/3300030246.fna.gz "tests/metagenome/APH(3')-IIIa.fa" -c 1000000 -p 8
```

### Matrix commands

```bash
pgr mat to-phylip tests/mat/IBPA.fa.tsv

pgr mat to-pair tests/mat/IBPA.phy

pgr mat format tests/mat/IBPA.phy

pgr mat subset tests/mat/IBPA.phy tests/mat/IBPA.list

hnsm dist seq tests/mat/IBPA.fa -k 7 -w 1 |
    pgr mat phylip stdin -o tests/mat/IBPA.71.phy

pgr mat compare tests/mat/IBPA.phy tests/mat/IBPA.71.phy --method all
```

---

## External Dependencies (Details)

Some subcommands in `pgr` invoke external executables. In particular:

- `pgr pl ucsc` depends on the UCSC kent-tools suite, including:
  - `faToTwoBit`, `axtChain`, `chainAntiRepeat`, `chainMergeSort`, `chainPreNet`,
    `chainNet`, `netSyntenic`, `netChainSubset`, `chainStitchId`, `netSplit`,
    `netToAxt`, `axtSort`, `axtToMaf`, `netFilter`, `netClass`, `chainSplit`.
- `pgr pl trf` depends on `trf` and `spanr`.
- `pgr pl rept` and `pgr pl ir` depend on `FastK`, `Profex`, and `spanr`.
- `pgr pl p2m` depends on `spanr`.
- `pgr fas refine` depends on a multiple sequence alignment tool:
  - `clustalw` (default), or `muscle`, or `mafft`.

Ensure these tools are installed and available in your `PATH` before running the corresponding pipelines.

---

## Help Text Style Guide

For contributors adding new subcommands, `pgr` uses a consistent style for CLI help text.

- **`about`**: Third-person singular (e.g., "Counts...", "Calculates...").
- **`after_help`**: Uses a raw string like `r###"..."###`.
  - **Description**: Detailed explanation.
  - **Notes**: Bullet points starting with `*`.
    - Standard note for `fa`/`fas`: `* Supports both plain text and gzipped (.gz) files`
    - Standard note for `fa`/`fas`: `* Reads from stdin if input file is 'stdin'`
    - Standard note for `twobit`: `* 2bit files are binary and require random access (seeking)`
    - Standard note for `twobit`: `* Does not support stdin or gzipped inputs`
  - **Examples**: Numbered list (`1.`, `2.`) with code blocks indented by 3 spaces.
- **Arguments**:
  - **Input**: `infiles` (multiple) or `infile` (single).
    - Help: `Input [FASTA|block FA|2bit] file(s) to process`.
  - **Output**: `outfile` (`-o`, `--outfile`).
    - Help: `Output filename. [stdout] for screen`.
- **Terminology**:
  - `pgr fa` → "FASTA"
  - `pgr fas` → "block FA"
  - `pgr twobit` → "2bit"

