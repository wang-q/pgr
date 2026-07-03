use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("window")
        .about("Splits sequences into overlapping windows")
        .after_help(
            r###"
This command splits sequences in a FASTA file into overlapping windows.

Header format:
    >seq_name:start-end

Notes:
* Coordinates are 1-based, inclusive.
* Windows containing only Ns are skipped.
* Output sequences are unwrapped (single line).

Coverage & Overlap:
* Theoretical Coverage = Window Length / Step Size.
* Examples:
  - --window-length 200 --step 100: 2x coverage (50% overlap).
  - --window-length 200 --step 200: 1x coverage (no overlap).
  - --window-length 200 --step 10:  20x coverage (95% overlap).

Splitting & Shuffling:
* --chunk N: Splits output into files with N records each (e.g., output.001.fa).
* --shuffle: Randomizes output records.
  - With --chunk: Buffers N records, shuffles, writes to file, clears buffer (Low memory).
  - Without --chunk: Buffers ALL records, shuffles, writes to single file (High memory).
* --chunk cannot be used with stdout.

Examples:
1. Split into 200bp windows with 100bp step:
   pgr fa window input.fa --window-length 200 --step 100

2. Split large file into chunks of 1M records with shuffling:
   pgr fa window input.fa --chunk 1000000 --shuffle -o split.fa

3. Use default settings (200bp window, 100bp step):
   pgr fa window input.fa

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input FASTA file to process",
        ))
        .arg(crate::cmd_pgr::args::window_arg_with_default(
            "200",
            "Window length",
        ))
        .arg(
            Arg::new("step")
                .long("step")
                .value_parser(value_parser!(usize))
                .default_value("100")
                .help("Step size"),
        )
        .arg(
            Arg::new("shuffle")
                .long("shuffle")
                .action(ArgAction::SetTrue)
                .help("Shuffle the output records (uses more memory)"),
        )
        .arg(crate::cmd_pgr::args::seed_arg(
            Some("42"),
            None,
            "Random seed for shuffling",
        ))
        .arg(
            Arg::new("chunk_records")
                .long("chunk-records")
                .value_parser(value_parser!(usize))
                .help("Split output into chunks of N records"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let len = *args.get_one::<usize>("window").unwrap();
    let step = *args.get_one::<usize>("step").unwrap();
    let shuffle = args.get_flag("shuffle");
    let seed = *args.get_one::<u64>("seed").unwrap();
    let chunk_size = args.get_one::<usize>("chunk_records").copied();
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    pgr::libs::fmt::fa::run_window(infile, len, step, shuffle, seed, chunk_size, outfile)
}
