use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for to-fas.
pub fn make_subcommand() -> Command {
    crate::cmd_pgr::args::add_poa_args(
        crate::cmd_pgr::args::add_query_args(crate::cmd_pgr::args::add_msa_flag(
            crate::cmd_pgr::args::add_fasta_tsv_arg(Command::new("to-fas")),
        )),
        false,
    )
    .arg(crate::cmd_pgr::args::outfile_arg())
    .about("Queries PAF index and outputs pairwise or multi-way block FASTA")
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

5. Write output to a file:
   pgr paf to-fas alignments.paf chr1:1000-5000 -f genomes.tsv -o out.fas

"###,
    )
}
/// Execute the to-fas command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let opts = crate::cmd_pgr::args::query_options_from_args(args);
    let (idx, all_results) = pgr::libs::paf::query::run_query(&opts)?;
    let fasta_tsv = args
        .get_one::<String>("fasta_tsv")
        .context("missing required argument: --fasta-tsv")?;
    let mut fasta_store = pgr::libs::paf::fasta::prepare_store(fasta_tsv, &idx)?;
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    if args.get_flag("msa") {
        let params = crate::cmd_pgr::args::get_poa_params(args);
        pgr::libs::paf::to_fas::write_msa_fas(
            &idx,
            &all_results,
            &mut fasta_store,
            params,
            &mut writer,
        )?;
    } else {
        pgr::libs::paf::to_fas::write_pairwise_fas(
            &idx,
            &all_results,
            &mut fasta_store,
            &mut writer,
        )?;
    }
    writer.flush()?;
    Ok(())
}
