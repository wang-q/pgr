use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};
use std::io::Write;
/// Build the clap subcommand for to-gfa.
pub fn make_subcommand() -> Command {
    crate::cmd_pgr::args::add_poa_args(
        crate::cmd_pgr::args::add_query_args(
            crate::cmd_pgr::args::add_fasta_tsv_arg(Command::new("to-gfa"))
                .arg(
                    Arg::new("crush")
                        .long("crush")
                        .action(ArgAction::SetTrue)
                        .help(
                            "Compress SNP bubbles (impg `crush` style; loses base-level ALT info)",
                        ),
                )
                .arg(crate::cmd_pgr::args::outfile_arg()),
        ),
        false,
    )
    .about("Queries PAF index and output local GFA via POA graph")
    .after_help(
        r###"
Queries a PAF file or saved index (same logic as `pgr paf query`) and
outputs a local GFA (v1.0) graph induced from a POA multiple sequence
alignment of each region.

For each region, all homologous fragments (target first, then each
query, '-' strand reverse-complemented) are fed into the POA engine.
The POA graph — nodes are bases, edges are adjacencies, paths trace
each input sequence — is then compacted and exported as GFA S/L/P.

Compaction (unchopping): linear stretches of single-base nodes with
no branching are merged into multi-base segments, reducing node count
by ~1 order of magnitude. SNP/indel bubbles remain as graph branches,
preserving the full variant structure.

--crush: optional impg `crush` style post-processing. Compresses simple
bubbles (nodes sharing the same in/out neighbors) into a single node,
keeping the highest-weight allele. Useful for SV overview graphs but
LOSES base-level ALT information — paths through ALTs are rewritten to
the kept allele. Off by default.

Each region produces an independent GFA block (node IDs restart at 1).
Multiple regions are separated by `# region: <name>` comment lines.

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
1. Single region to local GFA:
   pgr paf to-gfa alignments.paf chr1:1000-5000 -f genomes.tsv

2. Multi-way GFA with transitive BFS:
   pgr paf to-gfa alignments.paf chr1:1000-5000 -t -f genomes.tsv

3. SV overview graph with bubbles crushed:
   pgr paf to-gfa alignments.paf chr1:1000-5000 -t -f genomes.tsv --crush

4. Batch query from BED regions:
   pgr paf to-gfa alignments.paf.idx -b regions.bed -f genomes.tsv

"###,
    )
}
/// Execute the to-gfa command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let opts = crate::cmd_pgr::args::query_options_from_args(args);
    let (idx, all_results) = pgr::libs::paf::query::run_query(&opts)?;
    let mut fasta_store =
        pgr::libs::paf::fasta::prepare_store(args.get_one::<String>("fasta_tsv").unwrap(), &idx)?;

    let params = crate::cmd_pgr::args::get_poa_params(args);
    let crush = args.get_flag("crush");

    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    pgr::libs::paf::poa_compact::output_gfa(
        &mut writer,
        &idx,
        &all_results,
        &mut fasta_store,
        params,
        crush,
    )?;

    writer.flush()?;
    Ok(())
}
