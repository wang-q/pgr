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
                .default_value("100")
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
    let curdir = std::env::current_dir()?;
    let pgr = std::env::current_exe()?.display().to_string();
    let tempdir = tempfile::Builder::new().prefix("pgr_rept_").tempdir()?;
    let tempdir_str = tempdir.path().to_str().unwrap();

    run_cmd!(info "==> Paths")?;
    run_cmd!(info "    \"pgr\"     = ${pgr}")?;
    run_cmd!(info "    \"curdir\"  = ${curdir:?}")?;
    run_cmd!(info "    \"tempdir\" = ${tempdir_str}")?;

    run_cmd!(info "==> Absolute paths")?;
    let abs_infile = intspan::absolute_path(args.get_one::<String>("infile").unwrap())?
        .display()
        .to_string();
    let abs_outfile = crate::cmd_pgr::pl::common::abs_path_or_stdout(outfile)?;

    //----------------------------
    // Ops
    //----------------------------
    run_cmd!(info "==> Switch to tempdir")?;
    std::env::set_current_dir(tempdir_str)?;

    run_cmd!(info "==> FastK")?;
    run_cmd!(
        FastK -p -k${opt_kmer} -Ngenome ${abs_infile}
    )?;

    run_cmd!(info "==> Process each chromosome")?;
    run_cmd!(
        ${pgr} fa size ${abs_infile} -o chr.sizes
    )?;

    let chrs = crate::cmd_pgr::pl::common::read_chr_names("chr.sizes")?;

    let re_prof: regex::Regex = regex::Regex::new(
        r"(?xi)
            (?<start>\d+)       # start
            \s*-\s*             # spacer
            (?<end>\d+)         # end
            \s*                 # spacer
            \((?<depth>\d+)\)   # depth
            ",
    )?;

    let rg_files = crate::cmd_pgr::pl::common::run_profex_per_chr(&chrs, &re_prof, Some(2))?;

    run_cmd!(info "==> Outputs")?;
    run_cmd!(
        spanr cover $[rg_files] |
            spanr span --op fill -n ${opt_fk} stdin |
            spanr span --op excise -n ${opt_min} stdin |
            spanr span --op fill -n ${opt_ff} stdin -o ${abs_outfile}
    )?;

    //----------------------------
    // Done
    //----------------------------
    std::env::set_current_dir(&curdir)?;

    Ok(())
}
