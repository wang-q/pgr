use clap::{Arg, ArgAction, ArgMatches, Command};
use cmd_lib::run_cmd;

/// Build the clap subcommand for e-kmer.
pub fn make_subcommand() -> Command {
    Command::new("e-kmer")
        .about("Identifies repeats against an external library (k-mer)")
        .after_help(
            r###"
This command identifies repeats in a genome against an external repeat library
(Dfam, RepBase, or TnCentral) using k-mer analysis, mimicking the
functionality of `RepeatMasker`.

* <repeat> is path to the fasta file containing the repeat library.
* <infile> is path to fasta file, .fa.gz is supported. Cannot be stdin.

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
        .arg(crate::cmd_pgr::args::kmer_arg_with_default("17"))
        .arg(crate::cmd_pgr::args::fill_kmer_arg())
        .arg(crate::cmd_pgr::args::min_len_arg_with_default(
            "300",
            "Minimum length of repetitive fragments",
        ))
        .arg(crate::cmd_pgr::args::fill_fragment_arg())
        .arg(
            Arg::new("keep_index")
                .long("keep-index")
                .action(ArgAction::SetTrue)
                .help("Keep the built repeat table next to the library for reuse"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the e-kmer command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    crate::cmd_pgr::args::ensure_outfile_distinct(
        outfile,
        [
            args.get_one::<String>("repeat").unwrap().as_str(),
            args.get_one::<String>("infile").unwrap().as_str(),
        ],
    )?;

    let opt_kmer = *args.get_one::<usize>("kmer").unwrap();
    let opt_fk = *args.get_one::<usize>("fill_kmer").unwrap();
    let opt_min = *args.get_one::<usize>("min_len").unwrap();
    let opt_ff = *args.get_one::<usize>("fill_fragment").unwrap();
    let keep_index = args.get_flag("keep_index");

    let ctx = pgr::libs::pl::PipelineCtx::new("pgr_rept_e_")?;

    run_cmd!(info "==> Absolute paths")?;
    let abs_repeat = ctx.abs_path(args.get_one::<String>("repeat").unwrap())?;
    let abs_infile = ctx.abs_path(args.get_one::<String>("infile").unwrap())?;
    let abs_outfile = pgr::libs::pl::abs_path_or_stdout(outfile)?;

    let _cwd_guard = ctx.enter()?;

    let opts = pgr::libs::pl::RepeatOpts {
        pgr: ctx.pgr.clone(),
        abs_infile,
        abs_outfile,
        opt_kmer,
        opt_fk,
        opt_min,
        opt_ff,
        abs_repeat: Some(abs_repeat),
        keep_index,
        min_depth: None,
    };

    pgr::libs::pl::run_repeat_pipeline(&opts)?;

    // Done

    Ok(())
}
