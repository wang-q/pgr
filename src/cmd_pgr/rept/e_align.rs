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
* `--min-identity` uses the gap-compressed identity (matches + repeat
  matches over aligned bases, excluding insert bases), unlike `pgr sd`
  whose block identity includes inserts.

* All operations are running in a tempdir and no intermediate files are retained.

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
                .default_value("31")
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
                .help("K-mers occurring at least this often on either side are skipped as seeds"),
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
    crate::cmd_pgr::args::ensure_outfile_distinct(
        outfile,
        [
            args.get_one::<String>("repeat").unwrap().as_str(),
            args.get_one::<String>("infile").unwrap().as_str(),
        ],
    )?;
    let keep_index = args.get_flag("keep_index");

    let ctx = pgr::libs::pl::PipelineCtx::new("pgr_rept_ea_")?;

    run_cmd!(info "==> Absolute paths")?;
    let abs_repeat = ctx.abs_path(args.get_one::<String>("repeat").unwrap())?;
    let abs_infile = ctx.abs_path(args.get_one::<String>("infile").unwrap())?;
    let abs_outfile = pgr::libs::pl::abs_path_or_stdout(outfile)?;

    let min_identity = *args.get_one::<f64>("min_identity").unwrap();
    anyhow::ensure!(
        min_identity > 0.0 && min_identity <= 1.0,
        "--min-identity must be in (0, 1]: {}",
        min_identity
    );
    let kmer = *args.get_one::<usize>("kmer").unwrap();
    let smer = *args.get_one::<usize>("smer").unwrap();
    let window = *args.get_one::<usize>("window").unwrap();
    let parallel = *args.get_one::<usize>("parallel").unwrap();
    anyhow::ensure!(kmer > 0, "--kmer must be positive: {}", kmer);
    anyhow::ensure!(smer > 0, "--smer must be positive: {}", smer);
    anyhow::ensure!(window > 0, "--window must be positive: {}", window);
    anyhow::ensure!(parallel > 0, "--parallel must be positive: {}", parallel);

    let _cwd_guard = ctx.enter()?;

    let opts = pgr::libs::pl::AlignRepeatOpts {
        pgr: ctx.pgr.clone(),
        abs_repeat,
        abs_infile,
        abs_outfile,
        keep_index,
        kmer,
        smer,
        window,
        freq: *args.get_one::<usize>("freq").unwrap(),
        min_shared: *args.get_one::<usize>("min_shared").unwrap(),
        min_identity,
        min_len: *args.get_one::<usize>("min_len").unwrap(),
        fill_fragment: *args.get_one::<usize>("fill_fragment").unwrap(),
        parallel,
    };

    pgr::libs::pl::run_align_repeat_pipeline(&opts)?;

    Ok(())
}
