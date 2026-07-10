# Supported File Formats

pgr handles the following file formats. Each page documents the
structure, coordinate system, and pgr implementation details.

Formats with both a format reference (here in `formats/`) and a command
document (in `docs/`) are listed with both links.

| Format   | Short description                            | Format reference | Command document |
|----------|----------------------------------------------|------------------|------------------|
| 2bit     | UCSC binary genome sequence format           | [twobit.md](twobit.md) | [docs/twobit.md](../twobit.md) |
| AXT      | UCSC pairwise alignment (blastz output)      | [axt.md](axt.md) | [docs/axt.md](../axt.md) |
| Block FA | Multiple alignment blocks (pgr fas)          | — | [docs/fas.md](../fas.md) |
| Chain    | UCSC chained alignment blocks                | [chain.md](chain.md) | [docs/chain.md](../chain.md) |
| CIGAR    | Run-length alignment operations              | [cigar.md](cigar.md) | — |
| Distance | PHYLIP, Pairwise + matrix structures         | [distance.md](distance.md) | [docs/dist.md](../dist.md) |
| FASTA    | Nucleotide/protein sequences                 | — | [docs/fa.md](../fa.md) |
| FASTQ    | Sequences with quality scores                | — | [docs/fq.md](../fq.md) |
| GFF      | Generic Feature Format                       | — | [docs/gff.md](../gff.md) |
| LAV      | BLASTZ local alignment view                  | [lav.md](lav.md) | [docs/lav.md](../lav.md) |
| LOC      | FASTA random-access location index           | [loc.md](loc.md) | — |
| MAF      | Multiple Alignment Format                    | [maf.md](../maf.md) | [docs/maf.md](../maf.md) |
| Net      | UCSC hierarchical alignment net              | [net.md](net.md) | [docs/net.md](../net.md) |
| Newick   | Phylogenetic tree format                     | — | [docs/nwk.md](../nwk.md) |
| PAF      | Pairwise mApping Format (12-column TSV)      | — | [docs/paf.md](../paf.md) |
| PSL      | UCSC pairwise sequence alignment             | [psl.md](psl.md) | [docs/psl.md](../psl.md) |
| pbit     | pgr population 2bit + delta archive          | [pbit.md](../pbit.md) | [docs/pbit.md](../pbit.md) |

## Coordinate conventions

Most UCSC-derived formats use **0-based, half-open** intervals `[start, end)` in
pgr's internal representation, regardless of the file-format convention.
Exceptions are noted on each format page.
