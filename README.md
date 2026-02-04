# pgr - Practical Genome Refiner

[![Build](https://github.com/wang-q/pgr/actions/workflows/build.yml/badge.svg)](https://github.com/wang-q/pgr/actions)
[![codecov](https://codecov.io/gh/wang-q/pgr/branch/master/graph/badge.svg?token=8toyNHCsVU)](https://codecov.io/gh/wang-q/pgr)
[![license](https://img.shields.io/github/license/wang-q/pgr)](https://github.com//wang-q/pgr)

<!-- TOC -->
* [pgr - Population Genomes Refiner](#pgr---population-genomes-refiner)
  * [Install](#install)
  * [Synopsis](#synopsis)
    * [`pgr help`](#pgr-help)
  * [Examples](#examples)
    * [Genomes](#genomes)
    * [Block FA files](#block-fa-files)
  * [Author](#author)
  * [License](#license)
<!-- TOC -->

## Install

Current release: 0.1.0

```bash
cargo install --path . --force #--offline

# test
cargo test -- --test-threads=1

# build under WSL 2
mkdir -p /tmp/cargo
export CARGO_TARGET_DIR=/tmp/cargo
cargo build

```

## Synopsis

### `pgr help`

```text
`pgr` - Practical Genome Refiner

Usage: pgr [COMMAND]

Commands:
  ms-to-dna  Convert ms output haplotypes (0/1) to DNA sequences (FASTA)
  axt        Manipulate AXT alignment files
  chain      Manipulate Chain alignment files
  lav        Manipulate LAV alignment files
  maf        Manipulate MAF alignment files
  net        Manipulate Net alignment files
  psl        Manipulate PSL alignment files
  pl         Run integrated pipelines
  2bit       Manage 2bit files
  fa         Manipulate FASTA files
  fas        Manipulate block FA files
  fq         Manipulate FASTQ files
  help       Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version

Subcommand groups:

* Sequences:
    * 2bit - 2bit query and extraction
    * fa   - FASTA operations: info, records, transform, indexing
    * fas  - Block FA operations: info, subset, transform, file, variation
    * fq   - FASTQ interleaving and conversion

* Genome alignments:
    * chain - Chain operations: sort, filter, transform, to-net
    * net   - Net operations: info, subset, transform, convert
    * axt   - AXT sorting and conversion
    * lav   - Convert to PSL
    * maf   - Convert to Block FA
    * psl   - PSL statistics, manipulation, and conversion

* Pipelines:
    * pl - Pipeline tools: p2m, trf, ir, rept, ucsc

```

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

### Genomes

* genomes

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

* `hnsm distance`

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

* plot

```bash
FastGA -v -pafx tests/pgr/sakai.fa.gz tests/pgr/mg1655.fa.gz > tmp.paf
FastGA -v -psl tests/pgr/sakai.fa.gz tests/pgr/mg1655.fa.gz > tmp.psl

pgr chain -t="" tests/pgr/mg1655.fa.gz tests/pgr/sakai.fa.gz tmp.psl > tmp.chain.maf
pgr chain --syn -t="" tests/pgr/mg1655.fa.gz tests/pgr/sakai.fa.gz tmp.psl > tmp.syn.maf

lastz <(gzip -dcf tests/pgr/mg1655.fa.gz) <(gzip -dcf tests/pgr/sakai.fa.gz) |
    lavToPsl stdin stdout \
    > tmp.lastz.psl
pgr chain --syn -t="" tests/pgr/mg1655.fa.gz tests/pgr/sakai.fa.gz tmp.lastz.psl > tmp.lastz.maf

#wgatools maf2paf tmp.maf -o - |
#    sed 's/sakai\.fa\.//g' |
#    sed 's/mg1655\.fa\.//g' \
#    > tmp.paf
#PAFtoALN tmp.paf tests/pgr/sakai.fa.gz tests/pgr/mg1655.fa.gz
#ALNplot tmp -p -n0

wgatools dotplot -f paf tmp.paf > tmp.html
wgatools dotplot tmp.chain.maf > tmp.chain.html
wgatools dotplot tmp.syn.maf > tmp.syn.html
wgatools dotplot tmp.lastz.maf > tmp.lastz.html

#FastGA -v -1:tmp tests/pgr/sakai.fa.gz tests/pgr/mg1655.fa.gz
#ALNplot tmp -p -n0

```

| ![paf.png](images/paf.png) | ![chain.png](images/chain.png) |
|:--------------------------:|:------------------------------:|
|            paf             |             chain              |

| ![syn.png](images/syn.png) | ![lastz.png](images/lastz.png) |
|:--------------------------:|:------------------------------:|
|            syn             |             lastz              |

* repeats

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

* RepeatMasker

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

multiz M=10 tests/multiz/S288cvsRM11_1a.maf     tests/multiz/S288cvsSpar.maf     1 out1 out2

```

### Proteomes


* Hypervector

```bash
hnsm dist hv tests/clust/IBPA.fa
#tests/clust/IBPA.fa     tests/clust/IBPA.fa     776     776     776     776     0.0000  1.0000  1.0000
hnsm dist seq tests/clust/IBPA.fa --merge
#tests/clust/IBPA.fa     tests/clust/IBPA.fa     763     763     763     763     0.0000  1.0000  1.0000

hnsm dist hv tests/genome/mg1655.pro.fa.gz
#tests/genome/mg1655.pro.fa.gz    tests/genome/mg1655.pro.fa.gz    1240734 1240734 1240734 1240734 0.0000  1.0000  1.0000
hnsm dist seq tests/genome/mg1655.pro.fa.gz --merge
#tests/genome/mg1655.pro.fa.gz    tests/genome/mg1655.pro.fa.gz    1267403 1267403 1267403 1267403 0.0000  1.0000  1.0000

hnsm dist hv tests/genome/mg1655.pro.fa.gz tests/genome/pao1.pro.fa.gz -k 7 -w 1
#tests/genome/mg1655.pro.fa.gz    tests/genome/pao1.pro.fa.gz      1240734 1733273 81195   2892811 0.4154  0.0281  0.0654
hnsm dist seq tests/genome/mg1655.pro.fa.gz tests/genome/pao1.pro.fa.gz -k 7 -w 1 --merge
#tests/genome/mg1655.pro.fa.gz    tests/genome/pao1.pro.fa.gz      1267403 1770832 60605   2977630 0.4602  0.0204  0.0478

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
# Sequences in matrices: 10 and 10
# Common sequences: 10
# Method  Score
# pearson 0.935803
# spearman        0.919631
# mae     0.113433
# cosine  0.978731
# jaccard 0.759106
# euclid  1.229844

```


## Deps

pgr 子命令所依赖的外部可执行程序：

- pgr pl ucsc : 依赖 UCSC kent-tools 套件。
  - 包括: faToTwoBit , axtChain , chainAntiRepeat , chainMergeSort , chainPreNet , chainNet , netSyntenic , netChainSubset , chainStitchId , netSplit , netToAxt , axtSort , axtToMaf , netFilter , netClass , chainSplit 。
- pgr pl trf : 依赖 trf , spanr 。
- pgr pl rept / pgr pl ir : 依赖 FastK , Profex , spanr 。
- pgr pl p2m : 依赖 spanr 。
- pgr fas refine : 依赖 clustalw (默认), 或 muscle , mafft 。

## Help text style guide

* **`about`**: Third-person singular (e.g., "Counts...", "Calculates...").
* **`after_help`**: Uses raw string `r###"..."###`.
    * **Description**: Detailed explanation.
    * **Notes**: Bullet points starting with `*`.
        * Standard note for `fa`/`fas`: `* Supports both plain text and gzipped (.gz) files`
        * Standard note for `fa`/`fas`: `* Reads from stdin if input file is 'stdin'`
        * Standard note for `twobit`: `* 2bit files are binary and require random access (seeking)`
        * Standard note for `twobit`: `* Does not support stdin or gzipped inputs`
    * **Examples**: Numbered list (`1.`, `2.`) with code blocks indented by 3 spaces.
* **Arguments**:
    * **Input**: `infiles` (multiple) or `infile` (single).
        * Help: `Input [FASTA|block FA|2bit] file(s) to process`.
    * **Output**: `outfile` (`-o`, `--outfile`).
        * Help: `Output filename. [stdout] for screen`.
* **Terminology**:
    * `pgr fa` -> "FASTA"
    * `pgr fas` -> "block FA"
    * `pgr twobit` -> "2bit"

## Author

Qiang Wang <wang-q@outlook.com>

## License

MIT.

Copyright by Qiang Wang.

Written by Qiang Wang <wang-q@outlook.com>, 2024-
