//! Shared clap argument builders for subcommands.

use clap::{builder, Arg, ArgAction, ArgMatches, Command};

use pgr::libs::paf::query::QueryOptions;
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

/// `-o/--outdir` with a custom default value (not stdout).
pub fn outdir_arg_with_default(val: &'static str) -> Arg {
    Arg::new("outdir")
        .long("outdir")
        .short('o')
        .num_args(1)
        .default_value(val)
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

/// Required positional `infile` argument (caller must provide, may pass "stdin").
pub fn infile_arg_required() -> Arg {
    infile_arg_required_with_help("Input filename. [stdin] for standard input")
}

/// Required positional `infile` argument with a custom help text.
/// Index is auto-assigned by clap — do not add `.index(N)` to other positionals
/// unless this is the only positional or all positionals use explicit indices.
pub fn infile_arg_required_with_help(help: &'static str) -> Arg {
    Arg::new("infile").required(true).num_args(1).help(help)
}

/// Standard positional `infiles` argument (one or more, required) at index 1.
///
/// `label` is the format name used in the help text (e.g. `"FASTA"`,
/// `"block FA"`, `"2bit"`). Use inline definition with a different `.index(N)`
/// when another positional precedes `infiles`.
pub fn infiles_arg(label: &str) -> Arg {
    Arg::new("infiles")
        .required(true)
        .num_args(1..)
        .index(1)
        .help(format!("Input {label} file(s) to process"))
}

/// Positional `infiles` argument at a custom index (required, 1 or more files).
/// Use when another positional precedes `infiles`.
pub fn infiles_arg_at(label: &str, index: usize) -> Arg {
    Arg::new("infiles")
        .required(true)
        .num_args(1..)
        .index(index)
        .help(format!("Input {label} file(s) to process"))
}

/// Positional `infiles` argument at index 1 with a custom help text.
/// Use when the default "Input {label} file(s) to process" pattern doesn't fit.
pub fn infiles_arg_with_help(help: &'static str) -> Arg {
    Arg::new("infiles")
        .required(true)
        .num_args(1..)
        .index(1)
        .help(help)
}

/// Positional `infiles` argument at index 1 with custom num_args and help.
/// Use when the default `1..` range doesn't fit (e.g., `2..`, `1..=2`, `1..=4`).
pub fn infiles_arg_with_numargs(
    help: &'static str,
    num_args: impl clap::builder::IntoResettable<builder::ValueRange>,
) -> Arg {
    Arg::new("infiles")
        .required(true)
        .num_args(num_args)
        .index(1)
        .help(help)
}

/// Positional `target` genome file argument (required, index 1).
pub fn target_genome_arg(help: &'static str) -> Arg {
    Arg::new("target")
        .required(true)
        .index(1)
        .num_args(1)
        .help(help)
}

/// Positional `query` genome file argument (required, index 2).
pub fn query_genome_arg(help: &'static str) -> Arg {
    Arg::new("query")
        .required(true)
        .index(2)
        .num_args(1)
        .help(help)
}

/// Standard `-i/--invert` flag for `some`-style subcommands (invert selection).
pub fn invert_arg() -> Arg {
    invert_arg_with_help("Invert selection: output sequences NOT in the list")
}

/// `-i/--invert` flag with a custom help text.
pub fn invert_arg_with_help(help: &'static str) -> Arg {
    Arg::new("invert")
        .long("invert")
        .short('i')
        .action(ArgAction::SetTrue)
        .help(help)
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

/// `-p/--parallel` with a custom default value.
pub fn parallel_arg_with_default(default: &'static str) -> Arg {
    Arg::new("parallel")
        .long("parallel")
        .short('p')
        .num_args(1)
        .default_value(default)
        .value_parser(clap::value_parser!(usize))
        .help("Number of threads for parallel processing")
}

/// `--no-ns` flag (output size without Ns).
pub fn no_ns_arg() -> Arg {
    Arg::new("no_ns")
        .long("no-ns")
        .action(ArgAction::SetTrue)
        .help("Output size without Ns")
}

/// `-U/--upper` flag (convert sequences to uppercase).
pub fn upper_arg() -> Arg {
    Arg::new("upper")
        .long("upper")
        .short('U')
        .action(ArgAction::SetTrue)
        .help("Convert sequences to uppercase")
}

/// `-d/--dash` flag (remove dashes from sequences).
pub fn dash_arg() -> Arg {
    Arg::new("dash")
        .long("dash")
        .short('d')
        .action(ArgAction::SetTrue)
        .help("Remove dashes '-'")
}

/// `--t-name` argument with an optional default value.
pub fn t_name_arg(default: Option<&'static str>) -> Arg {
    let arg = Arg::new("t_name")
        .long("t-name")
        .num_args(1)
        .help("Custom name for the target genome");
    match default {
        Some(d) => arg.default_value(d),
        None => arg,
    }
}

/// `--q-name` argument with an optional default value.
pub fn q_name_arg(default: Option<&'static str>) -> Arg {
    let arg = Arg::new("q_name")
        .long("q-name")
        .num_args(1)
        .help("Custom name for the query genome");
    match default {
        Some(d) => arg.default_value(d),
        None => arg,
    }
}

/// `--seed` argument (u64) with an optional default, short flag, and help text.
pub fn seed_arg(default: Option<&'static str>, short: Option<char>, help: &'static str) -> Arg {
    let arg = Arg::new("seed")
        .long("seed")
        .num_args(1)
        .value_parser(clap::value_parser!(u64))
        .help(help);
    let arg = match default {
        Some(d) => arg.default_value(d),
        None => arg,
    };
    match short {
        Some(c) => arg.short(c),
        None => arg,
    }
}

/// `--name-prefix` argument with an optional default value.
pub fn name_prefix_arg(default: Option<&'static str>) -> Arg {
    let arg = Arg::new("name_prefix").long("name-prefix").num_args(1);
    match default {
        Some(d) => arg.default_value(d).help("Prefix of record names"),
        None => arg.help("Add prefix to sequence names"),
    }
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

/// Extract PAF query options from clap matches added by [`add_query_args`].
pub fn query_options_from_args(args: &ArgMatches) -> QueryOptions {
    QueryOptions {
        infile: args.get_one::<String>("infile").unwrap().clone(),
        region: args.get_one::<String>("region").cloned(),
        bed_regions: args.get_one::<String>("bed_regions").cloned(),
        transitive: args.get_flag("transitive"),
        max_depth: *args.get_one::<u16>("max_depth").unwrap(),
        min_len: *args.get_one::<usize>("min_len").unwrap(),
        min_dist: *args.get_one::<i32>("min_dist").unwrap(),
        min_identity: *args.get_one::<f64>("min_identity").unwrap(),
        min_output_len: *args.get_one::<i32>("min_output_len").unwrap(),
        merge_distance: *args.get_one::<i32>("merge_distance").unwrap(),
        min_degree: *args.get_one::<usize>("min_degree").unwrap(),
        min_chain_length: *args.get_one::<i32>("min_chain_length").unwrap(),
        subset_list: args.get_one::<String>("subset_list").cloned(),
        syntenic_filter: args.get_one::<String>("syntenic_filter").cloned(),
        fasta_tsv: args
            .try_get_one::<String>("fasta_tsv")
            .ok()
            .flatten()
            .cloned(),
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

/// Standard `--name` argument for fas subcommands (species name selector).
pub fn fas_name_arg(help: &'static str) -> Arg {
    Arg::new("name").long("name").num_args(1).help(help)
}

/// Standard `-g/--genome` argument for fas subcommands (reference genome FA file).
pub fn genome_arg() -> Arg {
    Arg::new("genome")
        .short('g')
        .long("genome")
        .required(true)
        .num_args(1)
        .help("Path to the reference genome FA file")
}

/// `-R/--required` argument (file with species names to keep).
pub fn required_species_list_arg() -> Arg {
    Arg::new("required")
        .long("required")
        .short('R')
        .required(true)
        .num_args(1)
        .help("File with a list of species names to keep, one per line")
}

/// `--suffix` argument with a custom default (file extension for output files).
pub fn suffix_arg(default: &'static str) -> Arg {
    Arg::new("suffix")
        .long("suffix")
        .num_args(1)
        .default_value(default)
        .help("File extension for output files")
}

/// `--engine` argument with parameterized possible values, default, and help.
/// Used by fas consensus (POA engine) and fas refine (MSA engine).
pub fn engine_arg(
    possible: &'static [&'static str],
    default: &'static str,
    help: &'static str,
) -> Arg {
    let values: Vec<builder::PossibleValue> = possible
        .iter()
        .map(|v| builder::PossibleValue::new(*v))
        .collect();
    Arg::new("engine")
        .long("engine")
        .num_args(1)
        .action(ArgAction::Set)
        .default_value(default)
        .value_parser(values)
        .help(help)
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

/// Standard `--regex` (`-x`) selector for nwk subcommands.
pub fn regex_arg() -> Arg {
    Arg::new("regex")
        .long("regex")
        .short('x')
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

/// `-M/--monophyly` flag with a custom help text.
pub fn monophyly_arg(help: &'static str) -> Arg {
    Arg::new("monophyly")
        .long("monophyly")
        .short('M')
        .action(ArgAction::SetTrue)
        .help(help)
}

/// Standard `-b/--bl` flag for nwk subcommands (keep branch lengths in output).
pub fn bl_arg() -> Arg {
    Arg::new("bl")
        .long("bl")
        .short('b')
        .action(ArgAction::SetTrue)
        .help("Keep branch lengths")
}

/// Standard `-l/--lca` argument for nwk subcommands (lowest common ancestor).
pub fn lca_arg() -> Arg {
    Arg::new("lca")
        .long("lca")
        .short('l')
        .num_args(1)
        .action(ArgAction::Append)
        .help("Lowest common ancestor of two nodes")
}

// ============================================================================
// clust subcommand builders
// ============================================================================

/// `--matrix` argument for clust commands (distance matrix file).
pub fn matrix_arg() -> Arg {
    Arg::new("matrix")
        .long("matrix")
        .num_args(1)
        .help("Distance matrix file")
}

/// Standard `--format` argument for clustering output.
pub fn format_arg() -> Arg {
    Arg::new("clust_format")
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

/// `--max-iter` argument (maximum iterations, default 100).
pub fn max_iter_arg() -> Arg {
    Arg::new("max_iter")
        .long("max-iter")
        .num_args(1)
        .default_value("100")
        .value_parser(clap::value_parser!(usize))
        .help("Maximum number of iterations")
}

/// `--method` argument for hierarchical clustering (default: ward).
pub fn clust_method_arg() -> Arg {
    Arg::new("clust_method")
        .long("method")
        .default_value("ward")
        .help("Clustering method (single, complete, average, weighted, centroid, median, ward)")
}

/// `--input-format` argument for clustering partition files (default: pair).
pub fn clust_input_format_arg() -> Arg {
    Arg::new("clust_input_format")
        .long("input-format")
        .value_parser([
            builder::PossibleValue::new("cluster"),
            builder::PossibleValue::new("pair"),
            builder::PossibleValue::new("long"),
        ])
        .default_value("pair")
        .help("Input format for partition files")
}

// ============================================================================
// mat subcommand builders
// ============================================================================

/// `--method` argument for matrix comparison (default: pearson).
/// Accepts comma-separated methods (e.g. "pearson,cosine") or "all".
/// Validation is done by the caller (each token checked against known methods).
pub fn mat_method_arg() -> Arg {
    Arg::new("mat_method")
        .long("method")
        .action(ArgAction::Set)
        .default_value("pearson")
        .help("Comparison method(s), comma-separated (all|pearson|spearman|mae|cosine|jaccard|euclid)")
}

/// `--format` argument for matrix output (default: full).
pub fn mat_format_arg() -> Arg {
    Arg::new("mat_format")
        .long("format")
        .action(ArgAction::Set)
        .value_parser([
            builder::PossibleValue::new("full"),
            builder::PossibleValue::new("lower"),
            builder::PossibleValue::new("strict"),
        ])
        .default_value("full")
        .help("Output format")
}

/// `--input-format` argument for matrix transform (default: phylip).
pub fn mat_input_format_arg() -> Arg {
    Arg::new("mat_input_format")
        .long("input-format")
        .default_value("phylip")
        .value_parser([
            builder::PossibleValue::new("phylip"),
            builder::PossibleValue::new("pair"),
        ])
        .help("Input format")
}

// ============================================================================
// paf subcommand builders (Command → Command transformers)
// ============================================================================

/// Add common query arguments to a clap Command.
/// Shared by `paf query`, `paf to-bed`, `paf to-gfa`, `paf to-vcf`,
/// `paf to-fas`, and `paf to-maf`.
pub fn add_query_args(cmd: Command) -> Command {
    cmd.arg(
        Arg::new("infile")
            .required(true)
            .index(1)
            .help("Input PAF file or .paf.idx index to query"),
    )
    .arg(
        Arg::new("region")
            .index(2)
            .help("Target region to query (e.g. chr1:1000-5000)"),
    )
    .arg(
        Arg::new("bed_regions")
            .long("bed-regions")
            .short('b')
            .num_args(1)
            .help("BED file with multiple regions for batch query (name start end per line)"),
    )
    .arg(
        Arg::new("transitive")
            .long("transitive")
            .short('t')
            .action(ArgAction::SetTrue)
            .help("Enable transitive BFS traversal"),
    )
    .arg(
        Arg::new("max_depth")
            .long("max-depth")
            .num_args(1)
            .default_value("2")
            .value_parser(clap::value_parser!(u16))
            .help("Maximum BFS depth (0 = unlimited, default: 2)"),
    )
    .arg(
        Arg::new("min_len")
            .long("min-len")
            .num_args(1)
            .default_value("10")
            .value_parser(clap::value_parser!(usize))
            .help("Minimum interval length to propagate (default: 10)"),
    )
    .arg(
        Arg::new("min_dist")
            .long("min-dist")
            .num_args(1)
            .default_value("10")
            .value_parser(clap::value_parser!(i32))
            .help("Minimum distance to merge adjacent intervals (default: 10)"),
    )
    .arg(
        Arg::new("min_identity")
            .long("min-identity")
            .num_args(1)
            .default_value("0.0")
            .value_parser(clap::value_parser!(f64))
            .help("Minimum gap-compressed identity (0.0-1.0, default: 0.0)"),
    )
    .arg(
        Arg::new("min_output_len")
            .long("min-output-len")
            .num_args(1)
            .default_value("0")
            .value_parser(clap::value_parser!(i32))
            .help("Minimum output interval length (default: 0 = no filter)"),
    )
    .arg(
        Arg::new("merge_distance")
            .long("merge-distance")
            .num_args(1)
            .default_value("0")
            .value_parser(clap::value_parser!(i32))
            .help("Merge adjacent output intervals within this distance (default: 0 = off)"),
    )
    .arg(
        Arg::new("min_degree")
            .long("min-degree")
            .num_args(1)
            .default_value("0")
            .value_parser(clap::value_parser!(usize))
            .help("Minimum distinct query sequences per region (default: 0 = off)"),
    )
    .arg(
        Arg::new("min_chain_length")
            .long("min-chain-length")
            .num_args(1)
            .default_value("0")
            .value_parser(clap::value_parser!(i32))
            .help("Minimum total aligned length per query (default: 0 = off)"),
    )
    .arg(
        Arg::new("subset_list")
            .long("subset-sequence-list")
            .num_args(1)
            .help("File with sequence names to include (one per line)"),
    )
    .arg(
        Arg::new("syntenic_filter")
            .long("syntenic-filter")
            .num_args(1)
            .help("UCSC chain file; drop query results whose query interval is not covered by any chain's query span (chain-level, both target and query name must match)"),
    )
}

/// Add the required `-f/--fasta-tsv` argument.
/// Shared by `paf to-gfa`, `paf to-vcf`, `paf to-fas`, and `paf to-maf`.
pub fn add_fasta_tsv_arg(cmd: Command) -> Command {
    cmd.arg(
        Arg::new("fasta_tsv")
            .long("fasta-tsv")
            .short('f')
            .required(true)
            .num_args(1)
            .help("TSV file: genome_name <tab> bgzf_fasta_path"),
    )
}

/// Add the optional `-f/--fasta-tsv` argument (for topology-only commands).
/// Shared by `paf stat` and `paf graph`.
pub fn add_optional_fasta_tsv_arg(cmd: Command) -> Command {
    cmd.arg(
        Arg::new("fasta_tsv")
            .long("fasta-tsv")
            .short('f')
            .num_args(1)
            .help("TSV file: genome_name <tab> bgzf_fasta_path (optional for topology-only mode)"),
    )
}

/// Add the `--min-var-len` argument (default 100).
/// Shared by `paf stat` and `paf graph`.
pub fn add_min_var_len_arg(cmd: Command) -> Command {
    cmd.arg(
        Arg::new("min_var_len")
            .long("min-var-len")
            .num_args(1)
            .default_value("100")
            .value_parser(clap::value_parser!(i32))
            .help("Minimum indel length to split at (default: 100)"),
    )
}

/// `--crush` flag for `paf to-gfa`.
pub fn crush_arg() -> Arg {
    Arg::new("crush")
        .long("crush")
        .action(ArgAction::SetTrue)
        .help("Compress SNP bubbles (impg 'crush' style; loses base-level ALT info)")
}

/// Add the `--msa` flag for POA-based multi-way output.
/// Shared by `paf to-fas` and `paf to-maf`.
pub fn add_msa_flag(cmd: Command) -> Command {
    cmd.arg(
        Arg::new("msa")
            .long("msa")
            .action(ArgAction::SetTrue)
            .help("Merge results per region into a multi-way block via POA"),
    )
}

// ============================================================================
// dist subcommand builders
// ============================================================================

/// `infiles` positional argument for dist subcommands (1 or 2 FA/list files).
pub fn pair_infiles_arg() -> Arg {
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
    kmer_arg_with_default("7")
}

/// `-k/--kmer` size argument with a custom default value.
pub fn kmer_arg_with_default(default: &'static str) -> Arg {
    Arg::new("kmer")
        .long("kmer")
        .short('k')
        .num_args(1)
        .default_value(default)
        .value_parser(clap::value_parser!(usize))
        .help("K-mer size")
}

/// `-w/--window` size argument (default: 1, for minimizers).
pub fn window_arg() -> Arg {
    window_arg_with_default("1", "Window size for minimizers")
}

/// `--sim` flag (convert distance to similarity).
pub fn sim_arg() -> Arg {
    Arg::new("sim")
        .long("sim")
        .action(ArgAction::SetTrue)
        .help("Convert distance to similarity (1 - distance)")
}

/// `--list-files` flag (treat infiles as list files).
pub fn list_arg() -> Arg {
    Arg::new("list_files")
        .long("list-files")
        .action(ArgAction::SetTrue)
        .help("Treat infiles as list files, where each line is a path to a sequence file")
}

/// Collect the `infiles` positional args as `&str` slices borrowing `args`.
pub fn collect_infiles(args: &ArgMatches) -> Vec<&str> {
    args.get_many::<String>("infiles")
        .unwrap()
        .map(|s| s.as_str())
        .collect::<Vec<_>>()
}

/// `--min-score` argument (f64) with a custom default value.
pub fn min_score_arg(default: &'static str) -> Arg {
    Arg::new("min_score")
        .long("min-score")
        .num_args(1)
        .default_value(default)
        .value_parser(clap::value_parser!(f64))
        .help("Minimum score threshold")
}

/// `--min-score` (f64) without default (optional threshold).
pub fn min_score_arg_optional(help: &'static str) -> Arg {
    Arg::new("min_score")
        .long("min-score")
        .num_args(1)
        .value_parser(clap::value_parser!(f64))
        .help(help)
}

/// `--min-len` (usize) without default.
pub fn min_len_arg() -> Arg {
    Arg::new("min_len")
        .long("min-len")
        .num_args(1)
        .value_parser(clap::value_parser!(usize))
        .help("Minimum length")
}

/// `--min-len` (usize) with a custom default and help text.
pub fn min_len_arg_with_default(default: &'static str, help: &'static str) -> Arg {
    Arg::new("min_len")
        .long("min-len")
        .num_args(1)
        .default_value(default)
        .value_parser(clap::value_parser!(usize))
        .help(help)
}

/// `--max-len` (usize) without default.
pub fn max_len_arg() -> Arg {
    Arg::new("max_len")
        .long("max-len")
        .num_args(1)
        .value_parser(clap::value_parser!(usize))
        .help("Maximum length")
}

/// `--max-score` (f64) without default (optional threshold).
pub fn max_score_arg_optional(help: &'static str) -> Arg {
    Arg::new("max_score")
        .long("max-score")
        .num_args(1)
        .value_parser(clap::value_parser!(f64))
        .help(help)
}

// ============================================================================
// Additional common builders
// ============================================================================

/// `-l/--line` sequence line length argument.
pub fn line_arg(default: Option<&'static str>) -> Arg {
    let arg = Arg::new("line")
        .long("line")
        .short('l')
        .num_args(1)
        .value_parser(clap::value_parser!(usize))
        .help("Sequence line length");
    match default {
        Some(d) => arg.default_value(d),
        None => arg,
    }
}

/// `-c/--chunk-size` argument (usize) with an optional default and custom help.
pub fn chunk_size_arg(default: Option<&'static str>, help: &'static str) -> Arg {
    let arg = Arg::new("chunk_size")
        .long("chunk-size")
        .short('c')
        .num_args(1)
        .value_parser(clap::value_parser!(usize))
        .help(help);
    match default {
        Some(d) => arg.default_value(d),
        None => arg,
    }
}

/// `-g/--gap` flag (only identify regions of N/n).
pub fn gap_arg() -> Arg {
    Arg::new("gap")
        .long("gap")
        .short('g')
        .action(ArgAction::SetTrue)
        .help("Only identify regions of N/n (gaps)")
}

/// `--no-mask` flag (do not apply sequence masking).
pub fn no_mask_arg() -> Arg {
    Arg::new("no_mask")
        .long("no-mask")
        .action(ArgAction::SetTrue)
        .help("Do not apply sequence masking")
}

/// Positional `ranges` argument (optional, index 2).
pub fn ranges_arg() -> Arg {
    Arg::new("ranges")
        .required(false)
        .index(2)
        .num_args(0..)
        .help("Ranges of interest")
}

/// `--replace-tsv` argument (required) for replace commands.
pub fn replace_tsv_arg() -> Arg {
    Arg::new("replace_tsv")
        .long("replace-tsv")
        .required(true)
        .num_args(1)
        .help("TSV file of original_name and replacement_name(s)")
}

/// `--runlist` argument (required) for region-based commands.
pub fn runlist_arg() -> Arg {
    Arg::new("runlist")
        .long("runlist")
        .required(true)
        .num_args(1)
        .help("JSON file of chromosome runlists")
}

/// `--mode` argument with possible values, a default, and a custom help text.
pub fn mode_arg(
    default: &'static str,
    possible: &'static [&'static str],
    help: &'static str,
) -> Arg {
    let values: Vec<builder::PossibleValue> = possible
        .iter()
        .map(|v| builder::PossibleValue::new(*v))
        .collect();
    Arg::new("mode")
        .long("mode")
        .num_args(1)
        .action(ArgAction::Set)
        .default_value(default)
        .value_parser(values)
        .help(help)
}

// ============================================================================
// fa subcommand additional builders
// ============================================================================

/// Positional `name_list` file argument for fa subcommands.
pub fn fa_name_list_arg(required: bool) -> Arg {
    Arg::new("name_list")
        .required(required)
        .index(2)
        .num_args(1)
        .help(if required {
            "File containing one sequence name per line"
        } else {
            "File containing one sequence name per line (optional)"
        })
}

// ============================================================================
// clust subcommand additional builders
// ============================================================================

/// `-k/--k` number of clusters argument.
pub fn k_arg() -> Arg {
    Arg::new("k")
        .long("k")
        .short('k')
        .num_args(1)
        .value_parser(clap::value_parser!(usize))
        .help("Number of clusters")
}

// ============================================================================
// chain subcommand builders
// ============================================================================

/// Positional `in_net` argument (required, auto-indexed).
pub fn in_net_arg() -> Arg {
    Arg::new("in_net")
        .required(true)
        .num_args(1)
        .help("Input net file")
}

/// Positional `in_chain` argument (required, auto-indexed).
pub fn in_chain_arg() -> Arg {
    Arg::new("in_chain")
        .required(true)
        .num_args(1)
        .help("Input chain file")
}

/// Positional `t_sizes` argument for chain subcommands (auto-indexed after infile).
pub fn chain_t_sizes_arg() -> Arg {
    Arg::new("t_sizes")
        .required(true)
        .num_args(1)
        .help("Target sizes file")
}

/// Positional `q_sizes` argument for chain subcommands (auto-indexed after t_sizes).
pub fn chain_q_sizes_arg() -> Arg {
    Arg::new("q_sizes")
        .required(true)
        .num_args(1)
        .help("Query sizes file")
}

/// `-t/--target-2bit` argument (target genome 2bit file).
pub fn target_2bit_arg() -> Arg {
    Arg::new("target_2bit")
        .long("target-2bit")
        .short('t')
        .required(true)
        .help("Target genome 2bit file")
}

/// `-q/--query-2bit` argument (query genome 2bit file).
pub fn query_2bit_arg() -> Arg {
    Arg::new("query_2bit")
        .long("query-2bit")
        .short('q')
        .required(true)
        .help("Query genome 2bit file")
}

/// Positional `psl` argument at index 3 (required).
pub fn psl_positional_arg(help: &'static str) -> Arg {
    Arg::new("psl")
        .required(true)
        .num_args(1)
        .index(3)
        .help(help)
}

/// `--incl-hap` flag (include haplotype sequences).
pub fn incl_hap_arg() -> Arg {
    Arg::new("incl_hap")
        .long("incl-hap")
        .action(ArgAction::SetTrue)
        .help("Include haplotype sequences")
}

/// `--gap-model` argument with parameterized default and possible values.
pub fn gap_model_arg(
    default: &'static str,
    possible: &'static [&'static str],
    help: &'static str,
) -> Arg {
    let values: Vec<builder::PossibleValue> = possible
        .iter()
        .map(|v| builder::PossibleValue::new(*v))
        .collect();
    Arg::new("gap_model")
        .long("gap-model")
        .num_args(1)
        .action(ArgAction::Set)
        .default_value(default)
        .value_parser(values)
        .help(help)
}

/// `--align-gap-open` (i32) argument (overrides --gap-model).
pub fn align_gap_open_arg() -> Arg {
    Arg::new("align_gap_open")
        .long("align-gap-open")
        .num_args(1)
        .value_parser(clap::value_parser!(i32))
        .allow_negative_numbers(true)
        .help("Alignment gap open cost (overrides --gap-model)")
}

/// `--align-gap-extend` (i32) argument (overrides --gap-model).
pub fn align_gap_extend_arg() -> Arg {
    Arg::new("align_gap_extend")
        .long("align-gap-extend")
        .num_args(1)
        .value_parser(clap::value_parser!(i32))
        .allow_negative_numbers(true)
        .help("Alignment gap extension cost (overrides --gap-model)")
}

/// `--score-scheme` argument (LASTZ format file or preset name like hoxd55).
pub fn score_scheme_arg() -> Arg {
    Arg::new("score_scheme")
        .long("score-scheme")
        .num_args(1)
        .help("Score scheme file (LASTZ format) or preset (e.g. hoxd55)")
}

// ============================================================================
// pl subcommand builders
// ============================================================================

/// `--fill-kmer` argument (default 2).
pub fn fill_kmer_arg() -> Arg {
    Arg::new("fill_kmer")
        .long("fill-kmer")
        .num_args(1)
        .default_value("2")
        .value_parser(clap::value_parser!(usize))
        .help("Fill holes between repetitive k-mers")
}

/// `--fill-fragment` argument (default 10).
pub fn fill_fragment_arg() -> Arg {
    Arg::new("fill_fragment")
        .long("fill-fragment")
        .num_args(1)
        .default_value("10")
        .value_parser(clap::value_parser!(usize))
        .help("Fill holes between repetitive fragments")
}

// ============================================================================
// Cross-domain shared builders (Round 4 additions)
// ============================================================================

/// `--color` argument (no short flag, optional default value).
pub fn color_arg(default: Option<&'static str>, help: &'static str) -> Arg {
    let arg = Arg::new("color").long("color").num_args(1);
    match default {
        Some(d) => arg.default_value(d),
        None => arg,
    }
    .help(help)
}

/// `--by-query` flag (split/sort on query instead of target).
pub fn by_query_arg(help: &'static str) -> Arg {
    Arg::new("by_query")
        .long("by-query")
        .action(ArgAction::SetTrue)
        .help(help)
}

/// `-C/--count` flag (count records/occurrences).
pub fn count_arg(help: &'static str) -> Arg {
    Arg::new("count")
        .long("count")
        .short('C')
        .action(ArgAction::SetTrue)
        .help(help)
}

/// `--syn` flag (synteny-related processing).
pub fn syn_arg(help: &'static str) -> Arg {
    Arg::new("syn")
        .long("syn")
        .action(ArgAction::SetTrue)
        .help(help)
}

/// `--type` argument for net subcommands (action varies: Set or Append).
pub fn net_type_arg(action: ArgAction, help: &'static str) -> Arg {
    Arg::new("type").long("type").action(action).help(help)
}

/// `-w/--window` size argument with a custom default and help text.
pub fn window_arg_with_default(default: &'static str, help: &'static str) -> Arg {
    Arg::new("window")
        .long("window")
        .short('w')
        .num_args(1)
        .default_value(default)
        .value_parser(clap::value_parser!(usize))
        .help(help)
}

// ============================================================================
// pbit subcommand builders
// ============================================================================

/// `--ref`/`-r` argument for `pbit create`.
pub fn pbit_ref_arg() -> Arg {
    Arg::new("ref")
        .long("ref")
        .short('r')
        .required(true)
        .num_args(1)
        .help("Reference FASTA file (plain or .gz)")
}

/// `-i/--infile` argument for `pbit create` / `pbit append`.
pub fn pbit_infiles_arg() -> Arg {
    Arg::new("infiles")
        .long("infile")
        .short('i')
        .required(false)
        .num_args(1)
        .action(ArgAction::Append)
        .help("Sample FASTA file(s) (plain or .gz)")
}

/// `--name` argument for `pbit create` / `pbit append` (TSV of sample info).
pub fn pbit_name_arg() -> Arg {
    Arg::new("name")
        .long("name")
        .num_args(1)
        .help("TSV file of `sample_name<TAB>fasta_path[<TAB>paf_path]` lines (overrides -i)")
}

/// `-p/--paf` argument for `pbit create` / `pbit append`.
pub fn pbit_paf_arg() -> Arg {
    Arg::new("paf")
        .long("paf")
        .short('p')
        .num_args(1)
        .action(ArgAction::Append)
        .help("PAF file(s) for CIGAR-driven encoding (paired with -i by order)")
}

/// `-s/--segment-size` argument for `pbit create`.
pub fn pbit_segment_size_arg() -> Arg {
    Arg::new("segment_size")
        .long("segment-size")
        .short('s')
        .num_args(1)
        .default_value("4096")
        .value_parser(clap::value_parser!(usize))
        .help("Reference segment size in bp (default: 4096)")
}

/// `-k/--kmer-len` argument for `pbit create`.
pub fn pbit_kmer_len_arg() -> Arg {
    Arg::new("kmer_len")
        .long("kmer-len")
        .short('k')
        .num_args(1)
        .default_value("15")
        .value_parser(clap::value_parser!(usize))
        .help("K-mer length for LZ-diff hashing (default: 15)")
}

/// `-l/--min-match-len` argument for `pbit create`.
pub fn pbit_min_match_len_arg() -> Arg {
    Arg::new("min_match_len")
        .long("min-match-len")
        .short('l')
        .num_args(1)
        .default_value("18")
        .value_parser(clap::value_parser!(u32))
        .help("Minimum match length for LZ-diff (default: 18)")
}

/// `--samples` flag for `pbit stat`.
pub fn pbit_samples_flag_arg() -> Arg {
    Arg::new("samples")
        .long("samples")
        .action(ArgAction::SetTrue)
        .help("List all sample names")
}

/// `--refs` flag for `pbit stat`.
pub fn pbit_refs_flag_arg() -> Arg {
    Arg::new("refs")
        .long("refs")
        .action(ArgAction::SetTrue)
        .help("List reference contigs (with segment counts)")
}

/// `--contigs` flag for `pbit stat`.
pub fn pbit_contigs_flag_arg() -> Arg {
    Arg::new("contigs")
        .long("contigs")
        .action(ArgAction::SetTrue)
        .help("List contigs per sample (or for a single sample with -s)")
}

/// `-s/--sample` argument for `pbit stat` / `pbit to-fa`.
pub fn pbit_sample_filter_arg(help: &'static str) -> Arg {
    Arg::new("sample")
        .long("sample")
        .short('s')
        .num_args(1)
        .help(help)
}
