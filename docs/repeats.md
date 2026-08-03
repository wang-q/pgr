# Repeat masking with pgr

Repeat masking finds repetitive regions of a genome and optionally hard/soft
masks them. This page covers repeat library downloads, running RepeatMasker
itself as a reference, and pgr's native repeat-detection commands.

## Repeat libraries

Repeat-masking workflows compare a genome against a library of known repeat
consensus sequences.

### TnCentral

TnCentral is a database of prokaryotic insertion sequences. Download the
complete IS library and prepare it for `pgr pl ir`:

```bash
# Download and unpack
curl -LO https://tncentral.ncc.unesp.br/api/download_blast/nc/tn_in_is
unzip -j tn_in_is 'tncentral_integrall_isfinder.fa'
gzip -9 -c 'tncentral_integrall_isfinder.fa' > tncentral.fa.gz

# Sanity check and quality filter
pgr fa size tests/pgr/tncentral.fa.gz
pgr dist seq tests/pgr/tncentral.fa.gz -k 17 -w 5 -p 8 |
    spanr filter stdin --ge 5:0.9
```

### RepBase

RepBase is the classic repeat database used by RepeatMasker. The distribution
tarball contains an EMBL-format library; convert it to FASTA with readseq:

```bash
curl -LO https://github.com/wang-q/ubuntu/releases/download/20190906/repeatmaskerlibraries-20140131.tar.gz
tar xvfz repeatmaskerlibraries-20140131.tar.gz Libraries/RepeatMaskerLib.embl

# https://sourceforge.net/projects/readseq/
java -jar ~/bin/readseq.jar -f fa Libraries/RepeatMaskerLib.embl
mv Libraries/RepeatMaskerLib.embl.fasta repbase.fa
gzip -9 -k repbase.fa
```

### Dfam

Dfam is a curated database of transposable element families. Its full
consensus FASTA can be used directly as the candidate set for masking
workflows (`pgr pl ir`, or a `pgr align pgi`-based masker):

*   Dfam website: <https://www.dfam.org/>
*   Dfam releases / downloads: <https://www.dfam.org/releases/>
*   The full consensus FASTA from the releases page works with `pgr pl ir`
    and with RepeatMasker's `-lib` option.

## RepeatMasker (reference)

RepeatMasker remains the reference annotation tool. Example run through a
singularity image, converting its `.out` to a GFF runlist for comparison:

```bash
singularity run ~/bin/repeatmasker_master.sif /app/RepeatMasker/RepeatMasker \
    ./genome.fa -xsmall -species "bacteria"

singularity run ~/bin/repeatmasker_master.sif /app/RepeatMasker/util/rmOutToGFF3.pl \
    ./genome.fa.out > mg1655.rm.gff

spanr gff tests/pgr/mg1655.rm.gff -o tests/pgr/mg1655.rm.json
```

## pgr native repeat detection

pgr detects repeats without RepeatMasker:

*   `pgr pl ir` — interspersed repeats, against a repeat library;
*   `pgr pl rept` — genome-internal repeats, no library needed;
*   `pgr pl trf` — tandem repeats, wrapping TRF.

All three emit runlist JSON ready for `pgr fa mask`:

```bash
pgr pl ir tests/pgr/tncentral.fa.gz tests/pgr/mg1655.fa.gz \
    > tests/pgr/mg1655.ir.json

spanr stat tests/pgr/mg1655.chr.sizes tests/pgr/mg1655.ir.json

pgr pl rept tests/pgr/mg1655.fa.gz \
    > tests/pgr/mg1655.rept.json

pgr pl trf tests/pgr/mg1655.fa.gz \
    > tests/pgr/mg1655.trf.json

spanr stat tests/pgr/mg1655.chr.sizes tests/pgr/mg1655.rm.json
spanr statop tests/pgr/mg1655.chr.sizes tests/pgr/mg1655.ir.json tests/pgr/mg1655.rm.json
```

## Notes

*   `pgr pl ir` accepts any repeat FASTA (Dfam, RepBase, TnCentral); see
    [pl.md](pl.md).
*   RepeatMasker normally uses a species-specific library via FamDB, or a
    custom library via `-lib`.
*   A native pgr masking plan (Dfam full library + `pgr align pgi` alignment
    -> runlist -> `fa mask`) is tracked in the design notes
    [notes/design/repeatmasker.md](../notes/design/repeatmasker.md).
