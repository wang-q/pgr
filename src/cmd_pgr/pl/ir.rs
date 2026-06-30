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
                .index(1)
                .help("The repeats database"),
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(2)
                .help("Input file to process"),
        )
        .arg(
            Arg::new("kmer")
                .long("kmer")
                .short('k')
                .num_args(1)
                .default_value("17")
                .value_parser(value_parser!(usize))
                .help("Size of the k-mer"),
        )
        .arg(
            Arg::new("fk")
                .long("fk")
                .num_args(1)
                .default_value("2")
                .value_parser(value_parser!(usize))
                .help("Fill holes between repetitive k-mers"),
        )
        .arg(
            Arg::new("min")
                .long("min")
                .num_args(1)
                .default_value("300")
                .value_parser(value_parser!(usize))
                .help("Minimum length of repetitive fragments"),
        )
        .arg(
            Arg::new("ff")
                .long("ff")
                .num_args(1)
                .default_value("10")
                .value_parser(value_parser!(usize))
                .help("Fill holes between repetitive fragments"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    let opt_kmer = *args.get_one::<usize>("kmer").unwrap();
    let opt_fk = *args.get_one::<usize>("fk").unwrap();
    let opt_min = *args.get_one::<usize>("min").unwrap();
    let opt_ff = *args.get_one::<usize>("ff").unwrap();

    //----------------------------
    // Paths
    //----------------------------
    let ctx = crate::cmd_pgr::pl::common::PipelineCtx::new("pgr_rm_")?;

    run_cmd!(info "==> Absolute paths")?;
    let abs_repeat = ctx.abs_path(args.get_one::<String>("repeat").unwrap())?;
    let abs_infile = ctx.abs_path(args.get_one::<String>("infile").unwrap())?;
    let abs_outfile = crate::cmd_pgr::pl::common::abs_path_or_stdout(outfile)?;

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

    let opts = crate::cmd_pgr::pl::common::RepeatOpts {
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

    crate::cmd_pgr::pl::common::run_repeat_pipeline(&opts)?;

    //----------------------------
    // Done
    //----------------------------
    ctx.leave()?;

    Ok(())
}
