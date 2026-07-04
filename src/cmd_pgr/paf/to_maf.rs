use clap::{ArgMatches, Command};
use std::io::Write;

use pgr::libs::paf::index::{PafIndex, QueryResult};
use pgr::libs::paf::msa_build::{build_msa_entries, build_pairwise_block, run_poa_msa};
use pgr::libs::poa::AlignmentParams;

// Output pairwise MAF blocks. Each QueryResult becomes one `a` block with two
// `s` lines (target first, query second). Sequences are fetched on demand via
// FastaStore; CIGAR is walked directly (no POA refinement).
#[allow(clippy::type_complexity)]
fn output_maf(
    idx: &PafIndex,
    all_results: &[((String, i32, i32), Vec<QueryResult>)],
    fasta_store: &mut pgr::libs::paf::fasta::FastaStore,
    writer: &mut dyn Write,
) -> anyhow::Result<()> {
    writeln!(writer, "##maf version=1")?;
    for (_, results) in all_results {
        for result in results {
            let blk = build_pairwise_block(idx, result, fasta_store)?;

            // size = number of non-gap bases
            let q_size = blk.q_aln.chars().filter(|c| *c != '-').count();
            let t_size = blk.t_aln.chars().filter(|c| *c != '-').count();

            writeln!(writer, "a")?;
            writeln!(
                writer,
                "s\t{0}\t{1}\t{2}\t+\t{3}\t{4}",
                blk.tname, blk.t_start, t_size, blk.t_src_size, blk.t_aln
            )?;
            writeln!(
                writer,
                "s\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}",
                blk.qname, blk.q_start_maf, q_size, blk.q_strand, blk.q_src_size, blk.q_aln
            )?;
            writeln!(writer)?;
        }
    }
    Ok(())
}

// Output multi-way MAF blocks via POA. For each region, collect target +
// all query sequences (queries RC'd if '-' strand), feed them into the POA
// engine, and emit one `a` block with N `s` lines. CIGAR is ignored.
#[allow(clippy::type_complexity)]
fn output_maf_msa(
    idx: &PafIndex,
    all_results: &[((String, i32, i32), Vec<QueryResult>)],
    fasta_store: &mut pgr::libs::paf::fasta::FastaStore,
    params: AlignmentParams,
    writer: &mut dyn Write,
) -> anyhow::Result<()> {
    writeln!(writer, "##maf version=1")?;
    for ((tname_region, _, _), results) in all_results {
        if results.is_empty() {
            continue;
        }

        let entries = build_msa_entries(idx, tname_region, results, fasta_store)?;
        let msa = run_poa_msa(&entries, params.clone());

        // Emit the MAF block.
        writeln!(writer, "a")?;
        for (e, aln) in entries.iter().zip(msa.iter()) {
            let size = aln.chars().filter(|c| *c != '-').count();
            writeln!(
                writer,
                "s\t{}\t{}\t{}\t{}\t{}\t{}",
                e.name, e.start, size, e.strand, e.src_size, aln
            )?;
        }
        writeln!(writer)?;
    }
    Ok(())
}
/// Build the clap subcommand for to-maf.
pub fn make_subcommand() -> Command {
    crate::cmd_pgr::args::add_poa_args(
        crate::cmd_pgr::args::add_query_args(crate::cmd_pgr::args::add_msa_flag(
            crate::cmd_pgr::args::add_fasta_tsv_arg(Command::new("to-maf")),
        )),
        false,
    )
    .arg(crate::cmd_pgr::args::outfile_arg())
    .about("Queries PAF index and output pairwise or multi-way MAF")
    .after_help(
        r###"
Queries a PAF file or saved index (same logic as `pgr paf query`) and
outputs MAF blocks.

Default mode (pairwise): each query result becomes one 2-sequence MAF
block restored directly from CIGAR. Alignments are assumed to be
already refined by chain/net — no POA refinement is performed.

--msa mode (multi-way): merge all query results of each region into a
single multi-sequence MAF block via POA. Sequences (target first, then
each query, '-' strand reverse-complemented) are fed into the POA
engine; CIGAR is ignored. Best used with --transitive to gather all
homologous fragments of a region.

-f/--fasta-tsv (required):
  TSV with two columns: genome_name <tab> bgzf_fasta_path
  Each genome_name must match a query/target name in the PAF index.
  A FASTA file may be referenced by multiple genome_names (multi-chrom).
  All genome names in the PAF index must be present in the TSV (strict
  validation — missing entries cause a hard error).

Notes:
* Input PAF files should contain cg:Z: tags for accurate projection
* Supports both plain text and gzipped (.gz) files (including BGZF)
* Reads from stdin if input file is 'stdin'

Examples:
1. Single region to pairwise MAF:
   pgr paf to-maf alignments.paf chr1:1000-5000 -f genomes.tsv

2. Multi-way MSA with transitive BFS:
   pgr paf to-maf alignments.paf chr1:1000-5000 -t --msa -f genomes.tsv

3. Batch query from BED regions:
   pgr paf to-maf alignments.paf.idx -b regions.bed -f genomes.tsv

"###,
    )
}
/// Execute the to-maf command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let opts = crate::cmd_pgr::args::query_options_from_args(args);
    let (idx, all_results) = pgr::libs::paf::query::run_query(&opts)?;
    let mut fasta_store =
        pgr::libs::paf::fasta::prepare_store(args.get_one::<String>("fasta_tsv").unwrap(), &idx)?;
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;
    if args.get_flag("msa") {
        let params = crate::cmd_pgr::args::get_poa_params(args);
        output_maf_msa(&idx, &all_results, &mut fasta_store, params, &mut writer)?;
    } else {
        output_maf(&idx, &all_results, &mut fasta_store, &mut writer)?;
    }
    Ok(())
}
