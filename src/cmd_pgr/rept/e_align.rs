use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};
use cmd_lib::run_cmd;

/// Build the clap subcommand for e-align.
pub fn make_subcommand() -> Command {
    Command::new("e-align")
        .about("Identifies repeats against an external library (alignment)")
        .after_help(
            r###"
This command identifies repeats in a genome against an external repeat library
(Dfam, RepBase, or TnCentral) by alignment, mimicking the masking behavior of
`RepeatMasker` without its annotation post-processing.

* <repeat> is path to the fasta file containing the repeat library.
* <infile> is path to fasta file, .fa.gz is supported. Cannot be stdin.

* The library is aligned to the genome with `pgr align pgi` (reference =
  genome, query = library). Alignment blocks are filtered by identity and
  length, merged into intervals, and written as a runlist JSON ready for
  `pgr fa mask`.

* Soft-masked (lowercase) genomes are detected and warned about: lowercase
  repeat regions fragment the alignment and drastically underestimate
  coverage, so uppercase the genome first (`tr a-z A-Z`) if warned.

* All operations are running in a tempdir and no intermediate files are retained.

* External dependencies
    * spanr

"###,
        )
        .arg(
            Arg::new("repeat")
                .required(true)
                .num_args(1)
                .help("The repeats database"),
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input file to process",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("kmer")
                .long("kmer")
                .short('k')
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("40")
                .help("k-mer size for indexing"),
        )
        .arg(
            Arg::new("smer")
                .long("smer")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("8")
                .help("Syncmer s-mer length for indexing"),
        )
        .arg(
            Arg::new("window")
                .long("window")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("5")
                .help("Syncmer window for indexing"),
        )
        .arg(
            Arg::new("freq")
                .long("freq")
                .short('f')
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("100")
                .help("Maximum k-mer frequency on either side to keep as seed"),
        )
        .arg(
            Arg::new("min_span")
                .long("min-span")
                .short('c')
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("50")
                .help("Minimum per-axis seed span (bp) for a chain"),
        )
        .arg(
            Arg::new("max_gap")
                .long("max-gap")
                .short('s')
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("1000")
                .help("Maximum bp gap between consecutive seeds in a chain"),
        )
        .arg(
            Arg::new("band")
                .long("band")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("128")
                .help("Diagonal band half-width (bp) around the chain mean"),
        )
        .arg(
            Arg::new("merge_gap")
                .long("merge-gap")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("5000")
                .help("Maximum gap (bp) between adjacent colinear chains to merge"),
        )
        .arg(
            Arg::new("min_shared")
                .long("min-shared")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("16")
                .help("Minimum shared seed length (bp)"),
        )
        .arg(
            Arg::new("workflow")
                .long("workflow")
                .value_parser(["greedy", "tube"])
                .num_args(1)
                .default_value("greedy")
                .help("Chaining workflow"),
        )
        .arg(
            Arg::new("min_identity")
                .long("min-identity")
                .value_parser(value_parser!(f64))
                .num_args(1)
                .default_value("0.70")
                .help("Minimum alignment identity"),
        )
        .arg(
            Arg::new("min_len")
                .long("min-len")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("50")
                .help("Minimum length of repetitive fragments"),
        )
        .arg(
            Arg::new("fill_fragment")
                .long("fill-fragment")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("10")
                .help("Fill holes between repetitive fragments"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("8")
                .help("Number of threads for parallel processing"),
        )
        .arg(
            Arg::new("keep_index")
                .long("keep-index")
                .action(ArgAction::SetTrue)
                .help("Keep the built pgi indexes next to the inputs for reuse"),
        )
}

/// Execute the e-align command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let keep_index = args.get_flag("keep_index");

    let ctx = pgr::libs::pl::PipelineCtx::new("pgr_rept_ea_")?;

    run_cmd!(info "==> Absolute paths")?;
    let abs_repeat = ctx.abs_path(args.get_one::<String>("repeat").unwrap())?;
    let abs_infile = ctx.abs_path(args.get_one::<String>("infile").unwrap())?;
    let abs_outfile = pgr::libs::pl::abs_path_or_stdout(outfile)?;

    let _cwd_guard = ctx.enter()?;

    let opts = pgr::libs::pl::AlignRepeatOpts {
        pgr: ctx.pgr.clone(),
        abs_repeat,
        abs_infile,
        abs_outfile,
        keep_index,
        kmer: *args.get_one::<usize>("kmer").unwrap(),
        smer: *args.get_one::<usize>("smer").unwrap(),
        window: *args.get_one::<usize>("window").unwrap(),
        freq: *args.get_one::<usize>("freq").unwrap(),
        min_span: *args.get_one::<usize>("min_span").unwrap(),
        max_gap: *args.get_one::<usize>("max_gap").unwrap(),
        band: *args.get_one::<usize>("band").unwrap(),
        merge_gap: *args.get_one::<usize>("merge_gap").unwrap(),
        min_shared: *args.get_one::<usize>("min_shared").unwrap(),
        workflow: args.get_one::<String>("workflow").unwrap().clone(),
        min_identity: *args.get_one::<f64>("min_identity").unwrap(),
        min_len: *args.get_one::<usize>("min_len").unwrap(),
        fill_fragment: *args.get_one::<usize>("fill_fragment").unwrap(),
        parallel: *args.get_one::<usize>("parallel").unwrap(),
    };

    pgr::libs::pl::run_align_repeat_pipeline(&opts)?;

    Ok(())
}
