use clap::*;
use cmd_lib::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("ir")
        .about("Identify interspersed repeats in a genome")
        .after_help(
            r###"
This command identifies interspersed repeats in a genome, mimicking the functionality of `RepeatMasker`.

* <repeat> is path to the fasta file containing repeats from Dfam, RepBase, or TnCentral。
* <infile> is path to fasta file, .fa.gz is supported. Cannot be stdin.

* All operations are running in a tempdir and no intermediate files are retained.

* External dependencies
    * FastK / Profex / Fastrm
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
        .arg(crate::cmd_pgr::args::kmer_arg_with_default("17"))
        .arg(crate::cmd_pgr::args::fill_kmer_arg())
        .arg(crate::cmd_pgr::args::min_len_arg_with_default(
            "300",
            "Minimum length of repetitive fragments",
        ))
        .arg(crate::cmd_pgr::args::fill_fragment_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    let opt_kmer = *args.get_one::<usize>("kmer").unwrap();
    let opt_fk = *args.get_one::<usize>("fill_kmer").unwrap();
    let opt_min = *args.get_one::<usize>("min_len").unwrap();
    let opt_ff = *args.get_one::<usize>("fill_fragment").unwrap();

    //----------------------------
    // Paths
    //----------------------------
    let ctx = pgr::libs::pl::PipelineCtx::new("pgr_rm_")?;

    run_cmd!(info "==> Absolute paths")?;
    let abs_repeat = ctx.abs_path(args.get_one::<String>("repeat").unwrap())?;
    let abs_infile = ctx.abs_path(args.get_one::<String>("infile").unwrap())?;
    let abs_outfile = pgr::libs::pl::abs_path_or_stdout(outfile)?;

    //----------------------------
    // Ops
    //----------------------------
    ctx.enter()?;

    let re_prof: regex::Regex = regex::Regex::new(
        r"(?xi)
            (?<start>\d+)       # start
            \s*-\s*             # spacer
            (?<end>\d+)         # end
            ",
    )?;

    let opts = pgr::libs::pl::RepeatOpts {
        pgr: ctx.pgr.clone(),
        abs_infile,
        abs_outfile,
        opt_kmer,
        opt_fk,
        opt_min,
        opt_ff,
        abs_repeat: Some(abs_repeat),
        re_prof,
        min_depth: None,
    };

    pgr::libs::pl::run_repeat_pipeline(&opts)?;

    //----------------------------
    // Done
    //----------------------------
    ctx.leave()?;

    Ok(())
}
