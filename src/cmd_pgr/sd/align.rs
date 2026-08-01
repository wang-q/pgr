//! `pgr sd align` — chain/net refinement of SD hits into PAF.

use anyhow::Context;
use clap::{Arg, ArgMatches, Command};
use cmd_lib::run_cmd;
use pgr::libs::paf::record::write_paf_record;
use std::io::Write;

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

    let ctx = pgr::libs::pl::PipelineCtx::new("pgr_sd_align_")?;
    let pgr = ctx.pgr.clone();
    let abs_genome = ctx.abs_path(genome)?;
    let abs_psl = ctx.abs_path(psl)?;
    let abs_outfile = if outfile == "stdout" {
        outfile.to_string()
    } else {
        ctx.abs_path(outfile)?
    };
    let _cwd_guard = ctx.enter()?;

    // Chain/net refinement: no --syn so non-syntenic (rearranged) SDs survive.
    run_cmd!(${pgr} pl chainnet ${abs_genome} ${abs_genome} ${abs_psl} -o chainnet_out)?;

    // Merge every MAF block into one PAF.
    let mut writer = pgr::writer(&abs_outfile)
        .with_context(|| format!("failed to open writer for {}", abs_outfile))?;
    for maf in pgr::libs::io::list_files_ext("chainnet_out", "maf") {
        let mut reader =
            pgr::reader(&maf).with_context(|| format!("failed to open MAF {}", maf))?;
        loop {
            let block = match pgr::libs::fmt::maf::next_maf_block(&mut reader) {
                Ok(b) => b,
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(e.into()),
            };
            if let Some(rec) = pgr::libs::paf::maf_import::maf_block_to_paf(&block)? {
                write_paf_record(&mut writer, &rec)?;
            }
        }
    }
    writer.flush()?;
    Ok(())
}
