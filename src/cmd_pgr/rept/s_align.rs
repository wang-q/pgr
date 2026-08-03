use clap::{value_parser, Arg, ArgMatches, Command};
use cmd_lib::run_cmd;

/// Build the clap subcommand for s-align.
pub fn make_subcommand() -> Command {
    Command::new("s-align")
        .about("Identifies repetitive regions by self alignment")
        .after_help(
            r###"
This command identifies repetitive regions of a genome by self-alignment,
porting the Cactus-style pipeline of `scripts/pgr-repeat.sh`: the genome is
split into overlapping windows, aligned back to itself with `lastz`, lifted
to genomic coordinates, and regions with coverage above a threshold are
written as a runlist JSON ready for `pgr fa mask`.

* <infile> is path to fasta file, .fa.gz is supported. Cannot be stdin.

* No repeat library is needed (self-comparison only).

* Soft-masked (lowercase) genomes are detected and warned about: lowercase
  repeat regions are skipped by lastz and underestimate coverage, so
  uppercase the genome first (`tr a-z A-Z`) if warned.

* All operations are running in a tempdir and no intermediate files are retained.

* External dependencies
    * lastz

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input file to process",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("window")
                .long("window")
                .short('w')
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("200")
                .help("Overlapping window length (bp)"),
        )
        .arg(
            Arg::new("step")
                .long("step")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("100")
                .help("Window step size (bp); 100 with window 200 gives 2x coverage"),
        )
        .arg(
            Arg::new("chunk_records")
                .long("chunk-records")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("10000")
                .help("Split window output into chunks of N records"),
        )
        .arg(
            Arg::new("preset")
                .long("preset")
                .value_parser(pgr::libs::lastz::preset_names())
                .num_args(1)
                .default_value("set01")
                .help("lastz parameter set"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("4")
                .help("Number of threads for parallel processing"),
        )
        .arg(
            Arg::new("min_depth")
                .long("min-depth")
                .short('m')
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("4")
                .help("Minimum alignment depth for a region to be kept"),
        )
}

/// Execute the s-align command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    let ctx = pgr::libs::pl::PipelineCtx::new("pgr_rept_sa_")?;

    run_cmd!(info "==> Absolute paths")?;
    let abs_infile = ctx.abs_path(args.get_one::<String>("infile").unwrap())?;
    let abs_outfile = pgr::libs::pl::abs_path_or_stdout(outfile)?;

    let _cwd_guard = ctx.enter()?;

    let opts = pgr::libs::pl::SelfAlignOpts {
        pgr: ctx.pgr.clone(),
        abs_infile,
        abs_outfile,
        window: *args.get_one::<usize>("window").unwrap(),
        step: *args.get_one::<usize>("step").unwrap(),
        chunk_records: *args.get_one::<usize>("chunk_records").unwrap(),
        preset: args.get_one::<String>("preset").unwrap().clone(),
        parallel: *args.get_one::<usize>("parallel").unwrap(),
        min_depth: *args.get_one::<usize>("min_depth").unwrap(),
    };

    pgr::libs::pl::run_self_align_pipeline(&opts)?;

    Ok(())
}
