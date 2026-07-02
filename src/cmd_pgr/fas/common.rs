//! Shared helpers for `pgr fas` subcommands.
//!
//! Provides a pipeline runner that reads FasBlock records from input files
//! and dispatches per-block processing either single-threaded or via a
//! crossbeam parallel pipeline (1 reader → N workers → 1 writer).

use pgr::libs::fmt::fas::{next_fas_block, run_parallel, FasBlock};
use std::io::Write;

/// Process FasBlock files either single-threaded or in parallel.
///
/// For each block read from `infiles`, calls `proc_block` to produce a string
/// chunk, and writes all chunks to `outfile`. When `parallel > 1`, uses a
/// crossbeam pipeline with `parallel` worker threads (output order may differ
/// from input order).
pub fn run_pipeline<F>(
    args: &clap::ArgMatches,
    parallel: usize,
    proc_block: F,
) -> anyhow::Result<()>
where
    F: Fn(&FasBlock) -> anyhow::Result<String> + Sync,
{
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    if parallel <= 1 {
        for infile in args.get_many::<String>("infiles").unwrap() {
            let mut reader = pgr::reader(infile)?;
            while let Ok(block) = next_fas_block(&mut reader) {
                let out_string = proc_block(&block)?;
                writer.write_all(out_string.as_ref())?;
            }
        }
    } else {
        let infiles: Vec<String> = args
            .get_many::<String>("infiles")
            .unwrap()
            .cloned()
            .collect();
        run_parallel(&infiles, parallel, &mut writer, &proc_block)?;
    }

    writer.flush()?;
    Ok(())
}
