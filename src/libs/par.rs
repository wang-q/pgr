//! Parallel pipeline primitives shared by `pgr dist` subcommands.
//!
//! Provides a writer thread + rayon pool pair, list/path resolution, two-set
//! entry loading, and a generic parallel pairwise iteration helper. None of
//! these depend on clap; the cmd layer extracts positional args and passes
//! them in.

use rayon::prelude::*;
use std::io::Write;
use std::thread::JoinHandle;

/// Spawn a writer thread draining a channel and configure the global rayon
/// pool with `num_threads`. Returns the sender and the writer join handle.
pub fn spawn_writer_and_pool(
    outfile: &str,
    num_threads: usize,
) -> anyhow::Result<(crossbeam::channel::Sender<String>, JoinHandle<()>)> {
    let (sender, receiver) = crossbeam::channel::bounded::<String>(256);

    let output = outfile.to_string();
    let writer_thread = std::thread::spawn(move || {
        let mut writer = crate::writer(&output).unwrap();
        for result in receiver {
            writer.write_all(result.as_bytes()).unwrap();
        }
    });

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()?;

    Ok((sender, writer_thread))
}

/// Resolve an infile to a list of paths. If `is_list` is true, read the file
/// as a one-path-per-line list; otherwise treat `infile` itself as the path.
pub fn resolve_paths(infile: &str, is_list: bool) -> anyhow::Result<Vec<String>> {
    if is_list {
        crate::libs::io::read_names::<Vec<String>>(infile)
    } else {
        Ok(vec![infile.to_string()])
    }
}

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
        let paths = resolve_paths(infiles[0], is_list)?;
        let entries = load_fn(&paths)?;
        Ok((entries.clone(), entries))
    } else {
        let paths1 = resolve_paths(infiles[0], is_list)?;
        let paths2 = resolve_paths(infiles[1], is_list)?;
        let entries1 = load_fn(&paths1)?;
        let entries2 = load_fn(&paths2)?;
        Ok((entries1, entries2))
    }
}

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
