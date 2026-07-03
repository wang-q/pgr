use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("interleave")
        .visible_alias("il")
        .about("Interleaves paired-end sequences")
        .after_help(
            r###"
This command interleaves paired-end sequences from one or two files.

Input modes:
* Two files: Interleave R1 and R2 files
* One file: Generate dummy R2 sequences (N's)

Features:
* Supports both FA and FQ formats
* Automatic format detection
* Custom read name prefix
* Custom starting index

Notes:
* Cannot read from stdin
* For FQ output, quality scores are:
  - Preserved from input FQ
  - Set to '!' (ASCII 33) for input FA
* Paired files must have same number of reads

Examples:
1. Interleave two FQ files:
   pgr fq interleave R1.fq R2.fq -o out.fq

2. Generate dummy pairs:
   pgr fq interleave R1.fa --name-prefix sample --start-index 1

3. Convert to FQ:
   pgr fq interleave R1.fa R2.fa --fq -o out.fq

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg_with_numargs(
            "Input FASTQ file(s)",
            1..=2,
        ))
        .arg(
            Arg::new("fq")
                .long("fq")
                .action(ArgAction::SetTrue)
                .help("Write FQ"),
        )
        .arg(crate::cmd_pgr::args::name_prefix_arg(Some("read")))
        .arg(
            Arg::new("start_index")
                .long("start-index")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("0")
                .help("Starting index"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    let is_out_fq = args.get_flag("fq");
    let opt_prefix = args.get_one::<String>("name_prefix").unwrap();
    let opt_start = *args.get_one::<usize>("start_index").unwrap();

    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();

    let _final_idx =
        pgr::libs::fmt::fq::interleave(&mut writer, &infiles, opt_prefix, opt_start, is_out_fq)?;

    Ok(())
}
