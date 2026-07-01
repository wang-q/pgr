# Supported File Formats

pgr handles the following file formats.  Each page documents the
structure, coordinate system, and pgr implementation details.

| Format   | Short description                            | Full reference |
|----------|----------------------------------------------|----------------|
| 2bit     | UCSC binary genome sequence format           | [twobit.md](twobit.md) |
| AXT      | UCSC pairwise alignment (blastz output)      | [axt.md](axt.md) |
| Chain    | UCSC chained alignment blocks                | [chain.md](chain.md), [docs/chain.md](../chain.md) |
| CIGAR    | Run-length alignment operations              | [cigar.md](cigar.md) |
| Distance | PHYLIP, Pairwise + matrix structures         | [distance.md](distance.md) |
| LAV      | BLASTZ local alignment view                  | [lav.md](lav.md) |
| LOC      | FASTA random-access location index           | [loc.md](loc.md) |
| Net      | UCSC hierarchical alignment net              | [net.md](net.md), [docs/net.md](../net.md) |
| PAF      | Pairwise mApping Format (12-column TSV)      | [docs/paf.md](../paf.md) |
| PSL      | UCSC pairwise sequence alignment             | [psl.md](psl.md), [docs/psl.md](../psl.md) |

Additional formats with dedicated top-level documentation:

| Format   | Document |
|----------|----------|
| FASTA    | [docs/fa.md](../fa.md) |
| FASTQ    | [docs/fq.md](../fq.md) |
| Block FA | [docs/fas.md](../fas.md) |
| GFF      | [docs/gff.md](../gff.md) |
| MAF      | [docs/maf.md](../maf.md) |
| Newick   | [docs/nwk.md](../nwk.md) |

## Coordinate conventions

Most UCSC-derived formats use **0-based, half-open** intervals `[start, end)` in
pgr's internal representation, regardless of the file-format convention.
Exceptions are noted on each format page.
