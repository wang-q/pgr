//! Shared clap argument builders for subcommands.

use clap::Arg;

/// Standard `-o/--outfile` argument defaulting to stdout.
pub fn outfile_arg() -> Arg {
    Arg::new("outfile")
        .long("outfile")
        .short('o')
        .num_args(1)
        .default_value("stdout")
        .help("Output filename. [stdout] for screen")
}

/// Standard `-r/--rgfile` argument (file of regions, one per line).
pub fn rgfile_arg() -> Arg {
    Arg::new("rgfile")
        .long("rgfile")
        .short('r')
        .num_args(1)
        .help("File of regions, one per line")
}

/// Standard `-p/--parallel` argument (number of threads, usize, default 1).
pub fn parallel_arg() -> Arg {
    Arg::new("parallel")
        .long("parallel")
        .short('p')
        .num_args(1)
        .default_value("1")
        .value_parser(clap::value_parser!(usize))
        .help("Number of threads for parallel processing")
}
