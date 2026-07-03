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

/// Standard `-i/--invert` flag for `some`-style subcommands (invert selection).
pub fn invert_arg() -> Arg {
    Arg::new("invert")
        .long("invert")
        .short('i')
        .action(ArgAction::SetTrue)
        .help("Invert selection: output sequences NOT in the list")
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

/// Extract PAF query options from clap matches added by [`add_query_args`].
pub fn query_options_from_args(args: &ArgMatches) -> QueryOptions {
    QueryOptions {
        infile: args.get_one::<String>("infile").unwrap().clone(),
        region: args.get_one::<String>("region").cloned(),
        bed_regions: args.get_one::<String>("bed_regions").cloned(),
        transitive: args.get_flag("transitive"),
        max_depth: *args.get_one::<u16>("max_depth").unwrap(),
        min_len: *args.get_one::<i32>("min_len").unwrap(),
        min_dist: *args.get_one::<i32>("min_dist").unwrap(),
        min_identity: *args.get_one::<f64>("min_identity").unwrap(),
        min_output_len: *args.get_one::<i32>("min_output_len").unwrap(),
        merge_distance: *args.get_one::<i32>("merge_distance").unwrap(),
        min_degree: *args.get_one::<usize>("min_degree").unwrap(),
        min_chain_length: *args.get_one::<i32>("min_chain_length").unwrap(),
        subset_list: args.get_one::<String>("subset_list").cloned(),
        syntenic_filter: args.get_one::<String>("syntenic_filter").cloned(),
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

/// Standard `-b/--bl` flag for nwk subcommands (keep branch lengths in output).
pub fn bl_arg() -> Arg {
    Arg::new("bl")
        .long("bl")
        .short('b')
        .action(ArgAction::SetTrue)
        .help("Keep branch lengths")
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
            .num_args(0)
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
            .value_parser(clap::value_parser!(i32))
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

/// Add the `--msa` flag for POA-based multi-way output.
/// Shared by `paf to-fas` and `paf to-maf`.
pub fn add_msa_flag(cmd: Command) -> Command {
    cmd.arg(
        Arg::new("msa")
            .long("msa")
            .num_args(0)
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
    Arg::new("kmer")
        .long("kmer")
        .short('k')
        .num_args(1)
        .default_value("7")
        .value_parser(clap::value_parser!(usize))
        .help("K-mer size")
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
