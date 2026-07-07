use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;

use pgr::libs::paf::to_maf::{write_msa_maf, write_pairwise_maf};

/// Build the clap subcommand for to-maf.
pub fn make_subcommand() -> Command {
    crate::cmd_pgr::args::add_poa_args(
        crate::cmd_pgr::args::add_query_args(crate::cmd_pgr::args::add_msa_flag(
            crate::cmd_pgr::args::add_fasta_tsv_arg(Command::new("to-maf")),
        )),
        false,
    )
    .arg(crate::cmd_pgr::args::outfile_arg())
    .about("Queries PAF index and outputs pairwise or multi-way MAF")
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
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    if args.get_flag("msa") {
        let params = crate::cmd_pgr::args::get_poa_params(args);
        write_msa_maf(&idx, &all_results, &mut fasta_store, params, &mut writer)?;
    } else {
        write_pairwise_maf(&idx, &all_results, &mut fasta_store, &mut writer)?;
    }
    writer.flush()?;
    Ok(())
}
