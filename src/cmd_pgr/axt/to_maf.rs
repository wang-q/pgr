use clap::{Arg, ArgAction, ArgMatches, Command};
/// Build the clap subcommand for to-maf.
pub fn make_subcommand() -> Command {
    Command::new("to-maf")
        .about("Converts from axt to maf format")
        .after_help(
            r###"
Where tSizes and qSizes are files that contain the sizes of the target and query sequences.
Very often this will be a chrom.sizes file.

Examples:
1. Convert axt to maf:
   pgr axt to-maf in.axt -t t.sizes -q q.sizes -o out.maf

2. Split output by target name:
   pgr axt to-maf in.axt -t t.sizes -q q.sizes --t-split -o out_dir
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg().help("Input AXT file. [stdin] for standard input"))
        .arg(crate::cmd_pgr::args::t_sizes_arg().required(true))
        .arg(crate::cmd_pgr::args::q_sizes_arg().required(true))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("q_prefix")
                .long("q-prefix")
                .help("Add prefix to start of query sequence name in maf")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("t_prefix")
                .long("t-prefix")
                .help("Add prefix to start of target sequence name in maf")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("t_split")
                .long("t-split")
                .help("Create a separate maf file for each target sequence. Output is a dir.")
                .action(ArgAction::SetTrue),
        )
}
/// Execute the to-maf command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = crate::cmd_pgr::args::get_infile(args);
    let t_sizes_path = args.get_one::<String>("t_sizes").unwrap();
    let q_sizes_path = args.get_one::<String>("q_sizes").unwrap();
    let output = crate::cmd_pgr::args::get_outfile(args);
    let q_prefix = args
        .get_one::<String>("q_prefix")
        .map(|s| s.as_str())
        .unwrap_or("");
    let t_prefix = args
        .get_one::<String>("t_prefix")
        .map(|s| s.as_str())
        .unwrap_or("");
    let t_split = args.get_flag("t_split");

    let t_sizes = pgr::read_sizes::<usize>(t_sizes_path)?;
    let q_sizes = pgr::read_sizes::<usize>(q_sizes_path)?;

    pgr::libs::fmt::axt::axt_to_maf(
        input, output, &t_sizes, &q_sizes, t_prefix, q_prefix, t_split,
    )
}
