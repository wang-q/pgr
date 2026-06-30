//! Shared helpers for `pgr fas` subcommands.
//!
//! Provides a pipeline runner that reads FasBlock records from input files
//! and dispatches per-block processing either single-threaded or via a
//! crossbeam parallel pipeline (1 reader → N workers → 1 writer).

use clap::Arg;
use pgr::libs::fmt::fas::FasBlock;
use std::io::Write;

/// Standard `-r/--required` argument (file of required species names).
pub fn required_arg() -> Arg {
    Arg::new("required")
        .long("required")
        .short('r')
        .required(true)
        .num_args(1)
        .help("Required: File with a list of species names to keep, one per line")
}

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
    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap())?;

    if parallel <= 1 {
        for infile in args.get_many::<String>("infiles").unwrap() {
            let mut reader = pgr::reader(infile)?;
            while let Ok(block) = pgr::libs::fmt::fas::next_fas_block(&mut reader) {
                let out_string = proc_block(&block)?;
                writer.write_all(out_string.as_ref())?;
            }
        }
    } else {
        run_parallel(args, parallel, &mut writer, &proc_block)?;
    }

    writer.flush()?;
    Ok(())
}

// Crossbeam pipeline: 1 reader → N workers → 1 writer.
fn run_parallel<W, F>(
    args: &clap::ArgMatches,
    parallel: usize,
    writer: &mut W,
    proc_block: &F,
) -> anyhow::Result<()>
where
    W: Write,
    F: Fn(&FasBlock) -> anyhow::Result<String> + Sync,
{
    let (snd1, rcv1) = crossbeam::channel::bounded::<FasBlock>(10);
    let (snd2, rcv2) = crossbeam::channel::bounded::<String>(10);

    crossbeam::scope(|s| {
        // Reader thread.
        s.spawn(|_| {
            for infile in args.get_many::<String>("infiles").unwrap() {
                let mut reader = match pgr::reader(infile) {
                    Ok(r) => r,
                    Err(_) => break,
                };
                while let Ok(block) = pgr::libs::fmt::fas::next_fas_block(&mut reader) {
                    if snd1.send(block).is_err() {
                        break;
                    }
                }
            }
            drop(snd1);
        });

        // Worker threads.
        for _ in 0..parallel {
            let (sendr, recvr) = (snd2.clone(), rcv1.clone());
            s.spawn(move |_| {
                for block in recvr.iter() {
                    if let Ok(out_string) = proc_block(&block) {
                        if sendr.send(out_string).is_err() {
                            break;
                        }
                    }
                }
            });
        }
        drop(snd2);

        // Writer thread (runs on this thread).
        for out_string in rcv2.iter() {
            if writer.write_all(out_string.as_ref()).is_err() {
                break;
            }
        }
    })
    .map_err(|_| anyhow::anyhow!("parallel pipeline failed (worker panic)"))?;

    Ok(())
}
