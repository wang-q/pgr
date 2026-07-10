use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;
/// Build the clap subcommand for to-vcf.
pub fn make_subcommand() -> Command {
    crate::cmd_pgr::args::add_poa_args(
        crate::cmd_pgr::args::add_query_args(crate::cmd_pgr::args::add_fasta_tsv_arg(
            Command::new("to-vcf"),
        )),
        false,
    )
    .arg(crate::cmd_pgr::args::outfile_arg())
    .about("Queries PAF index and outputs multi-way VCF via POA MSA")
    .after_help(
        r###"
Queries a PAF file or saved index (same logic as `pgr paf query`) and
outputs a VCF file with substitutions and indels called from a POA
multiple sequence alignment.

For each region, all homologous fragments (target first, then each
query, '-' strand reverse-complemented) are fed into the POA engine to
produce a multi-way MSA. Three variant classes are emitted:

* SNP: single target non-gap column where >=1 query differs. REF is
  the target base, ALT are the distinct non-REF bases.
* INS: consecutive target gap columns. REF is the 1bp anchor (target
  base just before the gap), ALT is anchor + inserted bases per sample.
* DEL: consecutive target non-gap columns where >=1 query has a gap.
  REF is the target segment, ALT is the per-sample non-gap concatenation.

GT fields encode each sample's allele (0=REF, 1..=N=ALT index, '.'=gap
or non-ACGT). INS and DEL are left-aligned against the reference: the
anchor position shifts leftward while the reference base before the
anchor equals the last base of every non-empty indel sequence. Partial
deletions (a sample deletes only part of the DEL region) fall back to
a best-effort allele and are not fully left-aligned.

Recommended with --transitive to gather all homologous fragments of
each region.

-f/--fasta-tsv (required):
  TSV with two columns: genome_name <tab> bgzf_fasta_path
  Each genome_name must match a query/target name in the PAF index.
  All genome names in the PAF index must be present in the TSV.

Notes:
* Input PAF files should contain cg:Z: tags (used for query projection)
* Supports both plain text and gzipped (.gz) files (including BGZF)
* Reads from stdin if input file is 'stdin'

Examples:
1. Single region to VCF:
   pgr paf to-vcf alignments.paf chr1:1000-5000 -f genomes.tsv

2. Multi-way VCF with transitive BFS:
   pgr paf to-vcf alignments.paf chr1:1000-5000 -t -f genomes.tsv

3. Batch query from BED regions:
   pgr paf to-vcf alignments.paf.idx -b regions.bed -f genomes.tsv

"###,
    )
}
/// Execute the to-vcf command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let opts = crate::cmd_pgr::args::query_options_from_args(args);
    let (idx, all_results, fasta_store_opt) = pgr::libs::paf::query::run_query(&opts)?;
    let mut fasta_store = fasta_store_opt.context("missing required argument: --fasta-tsv")?;

    let params = crate::cmd_pgr::args::get_poa_params(args);

    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    pgr::libs::paf::vcf::output_vcf(&mut writer, &idx, &all_results, &mut fasta_store, params)?;

    writer.flush()?;
    Ok(())
}
