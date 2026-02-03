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
  pipeline  UCSC chain/net pipeline
  ir        Identify interspersed repeats in a genome
  rept      Identify repetitive regions in a genome
  trf       Identify tandem repeats in a genome
  ms2dna    Convert ms output haplotypes (0/1) to DNA sequences (FASTA)
  axt       Axt tools
  chain     Chain tools
  lav       LAV tools
  net       Net tools
  psl       Psl tools
  2bit      2bit tools
  fa        Fasta tools
  fq        Fastq tools
  help      Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version

Subcommand groups:

* Fasta files
    * info: size / count / masked / n50
    * records: one / some / order / split
    * transform: replace / rc / filter / dedup / mask / sixframe
    * indexing: gz / range / prefilter

* Genome alignments:
    * chain
    * net
    * axt
    * lav
    * psl
    * 2bit
    * fa
    * fq

* Repeats:
    * ir / rept / trf

```

## Examples

### 2bit

```bash
# pgr fa to-2bit tests/fasta/ufasta.fa -o tests/fasta/ufasta.2bit
faToTwoBit tests/genome/mg1655.fa.gz tests/genome/mg1655.2bit

pgr 2bit size tests/genome/mg1655.2bit

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

## Deps

pgr 子命令所依赖的外部可执行程序：

- pgr pl ucsc : 依赖 UCSC kent-tools 套件。
  - 包括: faToTwoBit , axtChain , chainAntiRepeat , chainMergeSort , chainPreNet , chainNet , netSyntenic , netChainSubset , chainStitchId , netSplit , netToAxt , axtSort , axtToMaf , netFilter , netClass , chainSplit 。
- pgr pl trf : 依赖 trf , spanr 。
- pgr pl rept / pgr pl ir : 依赖 FastK , Profex , spanr 。
- pgr pl p2m : 依赖 spanr 。
- pgr fas refine : 依赖 clustalw (默认), 或 muscle , mafft 。

## Author

Qiang Wang <wang-q@outlook.com>

## License

MIT.

Copyright by Qiang Wang.

Written by Qiang Wang <wang-q@outlook.com>, 2024-
