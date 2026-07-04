use clap::{ArgMatches, Command};

use pgr::libs::paf::index::PafIndex;
use pgr::libs::paf::msa_build::{build_msa_entries, build_pairwise_block, run_poa_msa};
use pgr::libs::paf::query::QueryGroup;
use pgr::libs::poa::AlignmentParams;

fn output_fas_pairwise(
    idx: &PafIndex,
    all_results: &[QueryGroup],
    fasta_store: &mut pgr::libs::paf::fasta::FastaStore,
) -> anyhow::Result<()> {
    for (_, results) in all_results {
        for result in results {
            let blk = build_pairwise_block(idx, result, fasta_store)?;
            println!(">{0}(+):{1}-{2}", blk.tname, blk.t_start + 1, blk.t_end);
            println!("{}", blk.t_aln);
            println!(
                ">{0}({1}):{2}-{3}",
                blk.qname,
                blk.q_strand,
                blk.q_start_fwd + 1,
                blk.q_end_fwd
            );
            println!("{}", blk.q_aln);
            println!();
        }
    }
    Ok(())
}

fn output_fas_msa(
    idx: &PafIndex,
    all_results: &[QueryGroup],
    fasta_store: &mut pgr::libs::paf::fasta::FastaStore,
    params: AlignmentParams,
) -> anyhow::Result<()> {
    for ((tname_region, _, _), results) in all_results {
        let entries = build_msa_entries(idx, tname_region, results, fasta_store)?;
        if entries.is_empty() {
            continue;
        }

        let msa = run_poa_msa(&entries, params.clone());

        for (e, aln) in entries.iter().zip(msa.iter()) {
            let size = aln.chars().filter(|c| *c != '-').count() as i32;
            println!(
                ">{0}({3}):{1}-{2}",
                e.name,
                e.start + 1,
                e.start + size,
                e.strand
            );
            println!("{}", aln);
        }
        println!();
    }
    Ok(())
}
/// Build the clap subcommand for to-fas.
pub fn make_subcommand() -> Command {
    crate::cmd_pgr::args::add_poa_args(
        crate::cmd_pgr::args::add_query_args(crate::cmd_pgr::args::add_msa_flag(
            crate::cmd_pgr::args::add_fasta_tsv_arg(Command::new("to-fas")),
        )),
        false,
    )
    .about("Queries PAF index and output pairwise or multi-way block FASTA")
    .after_help(
        r###"
Queries a PAF file or saved index (same logic as `pgr paf query`) and
outputs block FASTA records.

Default mode (pairwise): each query result becomes a block of two FASTA
records (target first, query second) restored directly from CIGAR.
Alignments are assumed to be already refined by chain/net — no POA
refinement is performed.

--msa mode (multi-way): merge all query results of each region into a
single multi-sequence block FASTA via POA. Sequences (target first, then
each query, '-' strand reverse-complemented) are fed into the POA
engine; CIGAR is ignored. Best used with --transitive to gather all
homologous fragments of a region.

Output format (per block):
  >seq_name(+):start-end
  ATGC--ATGC
  >seq_name(-):start-end
  ATGCAT--GC
  (blank line)

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
* Output is compatible with `pgr fas to-vcf`

Examples:
1. Single region to pairwise FAS:
   pgr paf to-fas alignments.paf chr1:1000-5000 -f genomes.tsv

2. Multi-way MSA with transitive BFS:
   pgr paf to-fas alignments.paf chr1:1000-5000 -t --msa -f genomes.tsv

3. Pipeline to VCF:
   pgr paf to-fas alignments.paf chr1:1000-5000 -t --msa -f genomes.tsv | pgr fas to-vcf

4. Batch query from BED regions:
   pgr paf to-fas alignments.paf.idx -b regions.bed -f genomes.tsv

"###,
    )
}
/// Execute the to-fas command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let opts = crate::cmd_pgr::args::query_options_from_args(args);
    let (idx, all_results) = pgr::libs::paf::query::run_query(&opts)?;
    let mut fasta_store =
        pgr::libs::paf::fasta::prepare_store(args.get_one::<String>("fasta_tsv").unwrap(), &idx)?;
    if args.get_flag("msa") {
        let params = crate::cmd_pgr::args::get_poa_params(args);
        output_fas_msa(&idx, &all_results, &mut fasta_store, params)
    } else {
        output_fas_pairwise(&idx, &all_results, &mut fasta_store)
    }
}
