//! Shared clap argument builders for subcommands.

use clap::{Arg, ArgMatches, Command};

use pgr::libs::poa::AlignmentParams;

/// Standard `-o/--outfile` argument defaulting to stdout.
pub fn outfile_arg() -> Arg {
    Arg::new("outfile")
        .long("outfile")
        .short('o')
        .num_args(1)
        .default_value("stdout")
        .help("Output filename. [stdout] for screen")
}

/// Standard `-o/--outdir` argument defaulting to stdout.
pub fn outdir_arg() -> Arg {
    Arg::new("outdir")
        .long("outdir")
        .short('o')
        .num_args(1)
        .default_value("stdout")
        .help("Output directory. [stdout] for screen")
}

/// Standard positional `infile` argument defaulting to stdin.
pub fn infile_arg() -> Arg {
    Arg::new("infile")
        .num_args(1)
        .index(1)
        .default_value("stdin")
        .help("Input filename. [stdin] for standard input")
}

/// Standard positional `infiles` argument (one or more, required).
///
/// `label` is the format name used in the help text (e.g. `"FASTA"`,
/// `"block FA"`, `"2bit"`).
pub fn infiles_arg(label: &str) -> Arg {
    Arg::new("infiles")
        .required(true)
        .num_args(1..)
        .index(1)
        .help(format!("Input {label} file(s) to process"))
}

/// Standard `-r/--rgfile` argument (file of regions, one per line).
pub fn rgfile_arg() -> Arg {
    Arg::new("rgfile")
        .long("rgfile")
        .short('r')
        .num_args(1)
        .help("File of regions, one per line")
}

/// Standard `-t/--t-sizes` argument (target chromosome sizes file).
pub fn t_sizes_arg() -> Arg {
    Arg::new("t_sizes")
        .long("t-sizes")
        .short('t')
        .num_args(1)
        .help("Target sizes file")
}

/// Standard `-q/--q-sizes` argument (query chromosome sizes file).
pub fn q_sizes_arg() -> Arg {
    Arg::new("q_sizes")
        .long("q-sizes")
        .short('q')
        .num_args(1)
        .help("Query sizes file")
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

/// Extract the `outfile` value from `args` as `&str`.
pub fn get_outfile(args: &ArgMatches) -> &str {
    args.get_one::<String>("outfile").unwrap()
}

/// Extract the `infile` value from `args` as `&str`.
pub fn get_infile(args: &ArgMatches) -> &str {
    args.get_one::<String>("infile").unwrap()
}

/// Collect region strings from `ranges` (positional, optional) and `rgfile`
/// (`-r/--rgfile`) arguments. Returns the combined list.
pub fn collect_ranges(args: &ArgMatches) -> anyhow::Result<Vec<String>> {
    let mut ranges: Vec<String> = if args.contains_id("ranges") {
        args.get_many::<String>("ranges")
            .unwrap()
            .cloned()
            .collect()
    } else {
        Vec::new()
    };
    if args.contains_id("rgfile") {
        let mut rgs = pgr::libs::io::read_names_as_vec(args.get_one::<String>("rgfile").unwrap())?;
        ranges.append(&mut rgs);
    }
    Ok(ranges)
}

/// Add POA scoring arguments (`--match`, `--mismatch`, `--gap-open`,
/// `--gap-extend`) to `cmd`. When `with_shorts` is true, also registers the
/// `-m`/`-n`/`-g`/`-e` short flags (used by `fas consensus`; paf commands
/// pass false because `-m` collides with `--max-depth`).
pub fn add_poa_args(cmd: Command, with_shorts: bool) -> Command {
    let mut match_arg = Arg::new("match")
        .long("match")
        .num_args(1)
        .default_value("5")
        .value_parser(clap::value_parser!(i32))
        .allow_negative_numbers(true)
        .help("POA match score (default: 5)");
    if with_shorts {
        match_arg = match_arg.short('m');
    }

    let mut mismatch_arg = Arg::new("mismatch")
        .long("mismatch")
        .num_args(1)
        .default_value("-4")
        .value_parser(clap::value_parser!(i32))
        .allow_negative_numbers(true)
        .help("POA mismatch score (default: -4)");
    if with_shorts {
        mismatch_arg = mismatch_arg.short('n');
    }

    let mut gap_open_arg = Arg::new("gap_open")
        .long("gap-open")
        .num_args(1)
        .default_value("-8")
        .value_parser(clap::value_parser!(i32))
        .allow_negative_numbers(true)
        .help("POA gap open penalty (default: -8)");
    if with_shorts {
        gap_open_arg = gap_open_arg.short('g');
    }

    let mut gap_extend_arg = Arg::new("gap_extend")
        .long("gap-extend")
        .num_args(1)
        .default_value("-6")
        .value_parser(clap::value_parser!(i32))
        .allow_negative_numbers(true)
        .help("POA gap extend penalty (default: -6)");
    if with_shorts {
        gap_extend_arg = gap_extend_arg.short('e');
    }

    cmd.arg(match_arg)
        .arg(mismatch_arg)
        .arg(gap_open_arg)
        .arg(gap_extend_arg)
}

/// Extract POA scoring parameters from `ArgMatches` into `AlignmentParams`.
pub fn get_poa_params(args: &ArgMatches) -> AlignmentParams {
    AlignmentParams {
        match_score: *args.get_one::<i32>("match").unwrap(),
        mismatch_score: *args.get_one::<i32>("mismatch").unwrap(),
        gap_open: *args.get_one::<i32>("gap_open").unwrap(),
        gap_extend: *args.get_one::<i32>("gap_extend").unwrap(),
    }
}
