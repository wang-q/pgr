//! Shared helpers for `pgr dist` subcommands (seq / hv / vector).
//!
//! These subcommands share the same shape: a writer thread + rayon pool, a
//! pair of input file sets, and a parallel pairwise iteration that batches
//! output lines through a channel. This module exposes clap argument builders
//! plus a thin wrapper around the parallel primitives in `libs::par`.

use clap::{builder, Arg, ArgAction, ArgMatches};

// Re-export parallel primitives so callers can keep using `common::*`.
pub use pgr::libs::par::{load_entries, load_two_sets, par_run_pairs, spawn_writer_and_pool};

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
// Infile handling (clap-specific)
// ============================================================================

/// Collect the `infiles` positional args as `&str` slices borrowing `args`.
pub fn collect_infiles(args: &ArgMatches) -> Vec<&str> {
    args.get_many::<String>("infiles")
        .unwrap()
        .map(|s| s.as_str())
        .collect::<Vec<_>>()
}
