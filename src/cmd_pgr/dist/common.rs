//! Shared helpers for `pgr dist` subcommands (seq / hv / vector).
//!
//! These subcommands share the same shape: a writer thread + rayon pool, a
//! pair of input file sets, and a parallel pairwise iteration that batches
//! output lines through a channel. This module factors out that boilerplate.

use clap::{builder, Arg, ArgAction, ArgMatches};
use rayon::prelude::*;
use std::thread::JoinHandle;

// ============================================================================
// clap argument builders
// ============================================================================

/// `infiles` positional argument (1 or 2 FA/list files).
pub fn infiles_arg() -> Arg {
    Arg::new("infiles")
        .required(true)
        .num_args(1..=2)
        .index(1)
        .help("Input FA/list file(s). [stdin] for standard input")
}

/// `--hasher` selector (rapid / fx / murmur / mod).
pub fn hasher_arg() -> Arg {
    Arg::new("hasher")
        .long("hasher")
        .action(ArgAction::Set)
        .value_parser([
            builder::PossibleValue::new("rapid"),
            builder::PossibleValue::new("fx"),
            builder::PossibleValue::new("murmur"),
            builder::PossibleValue::new("mod"),
        ])
        .default_value("rapid")
        .help("Hash algorithm to use")
}

/// `-k/--kmer` size argument.
pub fn kmer_arg() -> Arg {
    Arg::new("kmer")
        .long("kmer")
        .short('k')
        .num_args(1)
        .default_value("7")
        .value_parser(clap::value_parser!(usize))
        .help("K-mer size")
}

/// `-w/--window` size argument.
pub fn window_arg() -> Arg {
    Arg::new("window")
        .long("window")
        .short('w')
        .num_args(1)
        .default_value("1")
        .value_parser(clap::value_parser!(usize))
        .help("Window size for minimizers")
}

/// `--sim` flag (convert distance to similarity).
pub fn sim_arg() -> Arg {
    Arg::new("sim")
        .long("sim")
        .action(ArgAction::SetTrue)
        .help("Convert distance to similarity (1 - distance)")
}

/// `--list` flag (treat infiles as list files).
pub fn list_arg() -> Arg {
    Arg::new("list")
        .long("list")
        .action(ArgAction::SetTrue)
        .help("Treat infiles as list files, where each line is a path to a sequence file")
}

// ============================================================================
// Writer thread + rayon thread pool
// ============================================================================

/// Spawn a writer thread draining a channel and configure the global rayon
/// pool with `num_threads`. Returns the sender and the writer join handle.
pub fn spawn_writer_and_pool(
    outfile: &str,
    num_threads: usize,
) -> anyhow::Result<(crossbeam::channel::Sender<String>, JoinHandle<()>)> {
    let (sender, receiver) = crossbeam::channel::bounded::<String>(256);

    let output = outfile.to_string();
    let writer_thread = std::thread::spawn(move || {
        let mut writer = pgr::writer(&output).unwrap();
        for result in receiver {
            writer.write_all(result.as_bytes()).unwrap();
        }
    });

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()?;

    Ok((sender, writer_thread))
}

// ============================================================================
// Infile handling
// ============================================================================

/// Collect the `infiles` positional args as `&str` slices borrowing `args`.
pub fn collect_infiles(args: &ArgMatches) -> Vec<&str> {
    args.get_many::<String>("infiles")
        .unwrap()
        .map(|s| s.as_str())
        .collect::<Vec<_>>()
}

/// Resolve an infile to a list of paths. If `is_list` is true, read the file
/// as a one-path-per-line list; otherwise treat `infile` itself as the path.
pub fn resolve_paths(infile: &str, is_list: bool) -> Vec<String> {
    if is_list {
        intspan::read_first_column(infile)
    } else {
        vec![infile.to_string()]
    }
}

// ============================================================================
// Entry loading
// ============================================================================

/// Load entries from a list of paths using a per-file loader.
pub fn load_entries<E, F>(paths: &[String], load_fn: F) -> anyhow::Result<Vec<E>>
where
    F: Fn(&str) -> anyhow::Result<Vec<E>>,
{
    let mut entries = Vec::new();
    for path in paths {
        let mut loaded = load_fn(path)?;
        entries.append(&mut loaded);
    }
    Ok(entries)
}

/// Load two entry sets for pairwise comparison.
///
/// With one infile: load it once and return `(entries.clone(), entries)` so
/// the caller can self-compare. With two infiles: load each independently.
/// `load_fn` receives the resolved path list for one set.
pub fn load_two_sets<E, F>(
    infiles: &[&str],
    is_list: bool,
    load_fn: F,
) -> anyhow::Result<(Vec<E>, Vec<E>)>
where
    E: Clone,
    F: Fn(&[String]) -> anyhow::Result<Vec<E>>,
{
    if infiles.len() == 1 {
        let paths = resolve_paths(infiles[0], is_list);
        let entries = load_fn(&paths)?;
        Ok((entries.clone(), entries))
    } else {
        let paths1 = resolve_paths(infiles[0], is_list);
        let paths2 = resolve_paths(infiles[1], is_list);
        let entries1 = load_fn(&paths1)?;
        let entries2 = load_fn(&paths2)?;
        Ok((entries1, entries2))
    }
}

// ============================================================================
// Parallel pairwise iteration
// ============================================================================

/// Iterate `entries1` x `entries2` in parallel (rayon), calling `pair_fn`
/// for each pair. If `pair_fn` returns `Some(line)`, the line is buffered
/// and flushed to `sender` every 1000 pairs (and at the end of each row).
pub fn par_run_pairs<E, F>(
    entries1: &[E],
    entries2: &[E],
    sender: &crossbeam::channel::Sender<String>,
    pair_fn: F,
) where
    E: Sync,
    F: Fn(&E, &E) -> Option<String> + Sync + Send,
{
    entries1.par_iter().for_each(|e1| {
        let mut lines = String::with_capacity(1024);
        for (i, e2) in entries2.iter().enumerate() {
            if let Some(out_string) = pair_fn(e1, e2) {
                lines.push_str(&out_string);
                if i > 1 && i % 1000 == 0 {
                    sender.send(lines.clone()).unwrap();
                    lines.clear();
                }
            }
        }
        if !lines.is_empty() {
            sender.send(lines).unwrap();
        }
    });
}
