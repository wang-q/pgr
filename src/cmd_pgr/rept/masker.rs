use clap::{value_parser, Arg, ArgMatches, Command};
use cmd_lib::run_cmd;

/// Build the clap subcommand for masker.
pub fn make_subcommand() -> Command {
    Command::new("masker")
        .about("Simulates RepeatMasker: rmblastn library search plus TRF simple repeats")
        .after_help(
            r###"
This command simulates RepeatMasker 4.2.4's `-lib` pipeline per batch:
TRF PERFECT (simple repeats, excised) -> rmblastn library search
(`general_search_parameters`) -> TRF DIVERGED (diverged simple repeats),
then writes a runlist JSON ready for `pgr fa mask`. RepeatMasker's annotation
post-processing (family/class, fragment re-joining, boundary refinement) is
not replicated.

* <repeat> is path to the fasta file containing the repeat library, .fa.gz
  is supported.
* <infile> is path to fasta file, .fa.gz is supported. Cannot be stdin.

* Search parameters (RepeatMasker 4.2.4 defaults): `-min_raw_gapped_score`
  = cutoff (225), `-word_size` by speed tier
  (slow/default/quick/rush = 8/9/11/13), `-gapopen 24 -gapextend 6`,
  `-mask_level 101`, `-complexity_adjust`, `-dust no`, xdrops 450/225/112,
  and the GC-keyed `20p##g.matrix` scoring matrix selected per fragment
  (RepeatMasker `chooseMatrices`).
* TRF stages use RepeatMasker's own parameters: PERFECT
  (2/7/7/80/10/50/10, copy > 4) then DIVERGED (2/3/5/75/20/33/7, copy > 5),
  with PERFECT simple repeats and library hits X-masked between stages like
  RepeatMasker's excise/mask flow.

* All operations run in a tempdir and no intermediate files are retained.

* External dependencies
    * makeblastdb and rmblastn (RMBlast >= 2.13), from $PATH or --rmblast-dir
    * trf, from $PATH

* Soft-masked (lowercase) genomes are warned about: rmblastn skips lowercase
  regions, so uppercase the genome first (`tr a-z A-Z`) if warned.

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
            Arg::new("cutoff")
                .long("cutoff")
                .value_parser(value_parser!(i32))
                .num_args(1)
                .default_value("225")
                .help("RepeatMasker cutoff score (default 225)"),
        )
        .arg(
            Arg::new("speed")
                .long("speed")
                .num_args(1)
                .default_value("default")
                .value_parser(["slow", "default", "quick", "rush"])
                .help("Search speed tier; sets -word_size (8/9/11/13)"),
        )
        .arg(
            Arg::new("matrix_gc")
                .long("matrix-gc")
                .value_parser(value_parser!(i64))
                .num_args(1)
                .help("Fixed GC percentage for scoring matrix selection (default: per chromosome)"),
        )
        .arg(
            Arg::new("min_len")
                .long("min-len")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("0")
                .help("Minimum length of repetitive fragments (0 = RepeatMasker raw hits)"),
        )
        .arg(
            Arg::new("fill_fragment")
                .long("fill-fragment")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("0")
                .help("Fill holes between repetitive fragments (0 = RepeatMasker raw hits)"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("8")
                .help("Total threads across rmblastn processes (4 per process, like RepeatMasker)"),
        )
        .arg(
            Arg::new("frag")
                .long("frag")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("60000")
                .help("Max fragment length before splitting (RepeatMasker -frag; 0 = whole chromosome)"),
        )
        .arg(
            Arg::new("rmblast_dir")
                .long("rmblast-dir")
                .num_args(1)
                .help(
                    "Directory containing makeblastdb and rmblastn (optional; falls back to $PATH)",
                ),
        )
}

/// Execute the masker command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    crate::cmd_pgr::args::ensure_outfile_distinct(
        outfile,
        [
            args.get_one::<String>("repeat").unwrap().as_str(),
            args.get_one::<String>("infile").unwrap().as_str(),
        ],
    )?;

    let cutoff = *args.get_one::<i32>("cutoff").unwrap();
    anyhow::ensure!(cutoff >= 0, "--cutoff must be non-negative: {}", cutoff);
    let speed = args.get_one::<String>("speed").unwrap().as_str();
    let word_size = match speed {
        "slow" => pgr::libs::rmblast::WORD_SIZE_TIERS[0],
        "default" => pgr::libs::rmblast::WORD_SIZE_TIERS[1],
        "quick" => pgr::libs::rmblast::WORD_SIZE_TIERS[2],
        "rush" => pgr::libs::rmblast::WORD_SIZE_TIERS[3],
        _ => unreachable!("speed validated by clap"),
    };
    let matrix_gc = args.get_one::<i64>("matrix_gc").copied();
    if let Some(gc) = matrix_gc {
        anyhow::ensure!(
            (0..=100).contains(&gc),
            "--matrix-gc must be in 0-100: {}",
            gc
        );
    }
    let parallel = *args.get_one::<usize>("parallel").unwrap();
    anyhow::ensure!(parallel > 0, "--parallel must be positive: {}", parallel);
    let frag = *args.get_one::<usize>("frag").unwrap();
    anyhow::ensure!(
        frag == 0 || frag >= 4000,
        "--frag must be 0 or >= 4000 (2x overlap): {}",
        frag
    );
    let rmblast_dir = args
        .get_one::<String>("rmblast_dir")
        .map(std::path::PathBuf::from);

    let ctx = pgr::libs::pl::PipelineCtx::new("pgr_rept_rm_")?;

    run_cmd!(info "==> Absolute paths")?;
    let abs_repeat = ctx.abs_path(args.get_one::<String>("repeat").unwrap())?;
    let abs_infile = ctx.abs_path(args.get_one::<String>("infile").unwrap())?;
    let abs_outfile = pgr::libs::pl::abs_path_or_stdout(outfile)?;

    let _cwd_guard = ctx.enter()?;

    let opts = pgr::libs::pl::MaskerOpts {
        pgr: ctx.pgr.clone(),
        abs_repeat,
        abs_infile,
        abs_outfile,
        cutoff,
        word_size,
        matrix_gc,
        min_len: *args.get_one::<usize>("min_len").unwrap(),
        fill_fragment: *args.get_one::<usize>("fill_fragment").unwrap(),
        parallel,
        frag,
        rmblast_dir,
    };

    pgr::libs::pl::run_masker_pipeline(&opts)?;

    Ok(())
}
