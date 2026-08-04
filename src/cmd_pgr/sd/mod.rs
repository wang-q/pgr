//! Subcommands for segmental duplication (SD) detection and analysis.

mod align;
mod cluster;
mod cover;
mod cross;
mod decompose;
mod run;
mod search;

use anyhow::Context;
use cmd_lib::run_cmd;
use pgr::libs::paf::record::write_paf_record;
use std::io::Write;

use clap::{ArgMatches, Command};

/// Build the `pgr sd` subcommand tree.
pub fn make_subcommand() -> Command {
    Command::new("sd")
        .about("Segmental duplication detection and analysis")
        .subcommand(align::make_subcommand())
        .subcommand(cluster::make_subcommand())
        .subcommand(cover::make_subcommand())
        .subcommand(cross::make_subcommand())
        .subcommand(decompose::make_subcommand())
        .subcommand(run::make_subcommand())
        .subcommand(search::make_subcommand())
}

/// Chain/net refine a PSL against target/query and merge the MAF into one PAF.
///
/// Runs `pgr pl chainnet` WITHOUT `--syn` (rearranged SDs survive) and merges
/// every output MAF block into a single PAF at `outfile` (stdout allowed).
pub(crate) fn chainnet_to_paf(
    target: &str,
    query: &str,
    psl: &str,
    outfile: &str,
) -> anyhow::Result<()> {
    let ctx = pgr::libs::pl::PipelineCtx::new("pgr_sd_refine_")?;
    let pgr = ctx.pgr.clone();
    let abs_target = ctx.abs_path(target)?;
    let abs_query = ctx.abs_path(query)?;
    let abs_psl = ctx.abs_path(psl)?;
    let abs_outfile = if outfile == "stdout" {
        outfile.to_string()
    } else {
        ctx.abs_path(outfile)?
    };
    let _cwd_guard = ctx.enter()?;

    run_cmd!(${pgr} pl chainnet ${abs_target} ${abs_query} ${abs_psl} -o chainnet_out)?;

    let mut writer = pgr::writer(&abs_outfile)
        .with_context(|| format!("failed to open writer for {}", abs_outfile))?;
    // Sort the per-contig MAF files so the merged PAF order is deterministic
    // across runs (read_dir order is filesystem-dependent).
    let mut maf_files = pgr::libs::io::list_files_ext("chainnet_out", "maf");
    maf_files.sort();
    for maf in maf_files {
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
            } else {
                log::warn!(
                    "skipping MAF block with {} component(s) (expected 2)",
                    block.components.len()
                );
            }
        }
    }
    writer.flush()?;
    Ok(())
}

/// Dispatch `pgr sd` subcommands.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("align", sub_matches)) => align::execute(sub_matches),
        Some(("cluster", sub_matches)) => cluster::execute(sub_matches),
        Some(("cover", sub_matches)) => cover::execute(sub_matches),
        Some(("cross", sub_matches)) => cross::execute(sub_matches),
        Some(("decompose", sub_matches)) => decompose::execute(sub_matches),
        Some(("run", sub_matches)) => run::execute(sub_matches),
        Some(("search", sub_matches)) => search::execute(sub_matches),
        _ => unreachable!("sd subcommand match"),
    }
}
