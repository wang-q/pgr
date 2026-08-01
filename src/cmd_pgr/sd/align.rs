//! `pgr sd align` — chain/net refinement of SD hits into PAF.

use clap::{Arg, ArgMatches, Command};

/// Build the clap subcommand for align.
pub fn make_subcommand() -> Command {
    Command::new("align")
        .about("Refines SD hits via chain/net and outputs PAF")
        .after_help(
            r###"
Runs the native chain-net-axt-maf pipeline (`pgr pl chainnet`, WITHOUT --syn so
rearranged SDs survive) on the putative hits from `pgr sd search`, then merges
the resulting MAF blocks into a single PAF for downstream cluster/decompose.

Notes:
* `target` and `query` are the same genome (self-alignment).
* Output PAF coordinates are 0-based half-open with cg:Z: CIGAR tags.

Examples:
1. Refine SD hits:
   pgr sd search genome.fa -o hits.psl
   pgr sd align genome.fa hits.psl -o hits.paf
"###,
        )
        .arg(
            Arg::new("genome")
                .index(1)
                .required(true)
                .help("Genome FASTA file (same as query and target)"),
        )
        .arg(
            Arg::new("psl")
                .index(2)
                .required(true)
                .help("Putative SD hits PSL file (from pgr sd search)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}
/// Execute the align command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let genome = args.get_one::<String>("genome").unwrap();
    let psl = args.get_one::<String>("psl").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    super::chainnet_to_paf(genome, genome, psl, outfile)
}
