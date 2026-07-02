use clap::*;
use cmd_lib::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("rept")
        .about("Identify repetitive regions in a genome")
        .after_help(
            r###"
This command identifies repetitive regions in a genome using k-mer analysis.

* <infile> is path to fasta file, .fa.gz is supported. Cannot be stdin.

* All operations are running in a tempdir and no intermediate files are retained.

* External dependencies
    * FastK / Profex / Fastrm
    * spanr

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Input file to process"),
        )
        .arg(
            Arg::new("kmer")
                .long("kmer")
                .short('k')
                .num_args(1)
                .default_value("17")
                .value_parser(value_parser!(usize))
                .help("K-mer size"),
        )
        .arg(
            Arg::new("fill_kmer")
                .long("fill-kmer")
                .num_args(1)
                .default_value("2")
                .value_parser(value_parser!(usize))
                .help("Fill holes between repetitive k-mers"),
        )
        .arg(
            Arg::new("min_len")
                .long("min-len")
                .num_args(1)
                .default_value("100")
                .value_parser(value_parser!(usize))
                .help("Minimum length of repetitive fragments"),
        )
        .arg(
            Arg::new("fill_fragment")
                .long("fill-fragment")
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
    let opt_fk = *args.get_one::<usize>("fill_kmer").unwrap();
    let opt_min = *args.get_one::<usize>("min_len").unwrap();
    let opt_ff = *args.get_one::<usize>("fill_fragment").unwrap();

    //----------------------------
    // Paths
    //----------------------------
    let ctx = crate::cmd_pgr::pl::common::PipelineCtx::new("pgr_rept_")?;

    run_cmd!(info "==> Absolute paths")?;
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
            \s*                 # spacer
            \((?<depth>\d+)\)   # depth
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
        abs_repeat: None,
        re_prof,
        min_depth: Some(2),
    };

    crate::cmd_pgr::pl::common::run_repeat_pipeline(&opts)?;

    //----------------------------
    // Done
    //----------------------------
    ctx.leave()?;

    Ok(())
}
