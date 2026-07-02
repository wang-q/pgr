//! Shared clap argument builders for subcommands.

use clap::{builder, Arg, ArgAction, ArgMatches, Command};

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

/// `-o/--outfile` with a custom default value.
pub fn outfile_arg_with_default(val: &'static str) -> Arg {
    Arg::new("outfile")
        .long("outfile")
        .short('o')
        .num_args(1)
        .default_value(val)
        .help("Output filename. [stdout] for screen")
}

/// `-o/--outfile` without default (optional). Caller must handle `None`.
pub fn outfile_arg_optional() -> Arg {
    Arg::new("outfile")
        .long("outfile")
        .short('o')
        .num_args(1)
        .help("Output filename. [stdout] for screen")
}

/// `-o/--outfile` required (no default).
pub fn outfile_arg_required() -> Arg {
    Arg::new("outfile")
        .long("outfile")
        .short('o')
        .num_args(1)
        .required(true)
        .help("Output filename")
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

/// `-o/--outdir` required (no default). For commands that must write to a directory.
pub fn outdir_arg_required() -> Arg {
    Arg::new("outdir")
        .long("outdir")
        .short('o')
        .num_args(1)
        .required(true)
        .help("Output directory")
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
        let mut rgs =
            pgr::libs::io::read_names::<Vec<String>>(args.get_one::<String>("rgfile").unwrap())?;
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

// ============================================================================
// fas subcommand builders
// ============================================================================

/// Standard `--outgroup` flag for fas subcommands.
pub fn outgroup_arg() -> Arg {
    Arg::new("outgroup")
        .long("outgroup")
        .action(ArgAction::SetTrue)
        .help("Indicates the presence of outgroups at the end of each block")
}

// ============================================================================
// nwk subcommand builders
// ============================================================================

/// Standard `--node` (`-n`) selector for nwk subcommands.
pub fn node_arg() -> Arg {
    Arg::new("node")
        .long("node")
        .short('n')
        .num_args(1)
        .action(ArgAction::Append)
        .help("Select nodes by exact name")
}

/// Standard `--name-list` (`-l`) selector for nwk subcommands.
pub fn name_list_arg() -> Arg {
    Arg::new("name_list")
        .long("name-list")
        .short('l')
        .num_args(1)
        .help("Select nodes from a name-list file")
}

/// Standard `--regex` (`-r`) selector for nwk subcommands.
pub fn regex_arg() -> Arg {
    Arg::new("regex")
        .long("regex")
        .short('r')
        .num_args(1)
        .action(ArgAction::Append)
        .help("Select nodes by regular expression (case insensitive)")
}

/// Standard `--descendants` (`-D`) flag for nwk subcommands.
pub fn descendants_arg() -> Arg {
    Arg::new("descendants")
        .long("descendants")
        .short('D')
        .action(ArgAction::SetTrue)
        .help("Include all descendants of selected internal nodes")
}

/// Standard `--internal` (`-I`) filter flag for nwk subcommands.
pub fn internal_arg() -> Arg {
    Arg::new("internal")
        .long("internal")
        .short('I')
        .action(ArgAction::SetTrue)
        .help("Don't print internal labels")
}

/// Standard `--leaf` (`-L`) filter flag for nwk subcommands.
pub fn leaf_arg() -> Arg {
    Arg::new("leaf")
        .long("leaf")
        .short('L')
        .action(ArgAction::SetTrue)
        .help("Don't print leaf labels")
}

// ============================================================================
// clust subcommand builders
// ============================================================================

/// Standard `--format` argument for clustering output.
pub fn format_arg() -> Arg {
    Arg::new("format")
        .long("format")
        .action(ArgAction::Set)
        .value_parser([
            builder::PossibleValue::new("cluster"),
            builder::PossibleValue::new("pair"),
        ])
        .default_value("cluster")
        .help("Output format for clustering results")
}

/// Standard `--same` argument. `default` varies by algorithm (mcl=1.0, dbscan/k-medoids=0.0).
pub fn same_arg(default: &'static str) -> Arg {
    Arg::new("same")
        .long("same")
        .num_args(1)
        .default_value(default)
        .value_parser(clap::value_parser!(f32))
        .help("Default score of identical element pairs")
}

/// Standard `--missing` argument. `default` varies by algorithm (mcl=0.0, dbscan/k-medoids=1.0).
pub fn missing_arg(default: &'static str) -> Arg {
    Arg::new("missing")
        .long("missing")
        .num_args(1)
        .default_value(default)
        .value_parser(clap::value_parser!(f32))
        .help("Default score of missing pairs")
}
