use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("count")
        .about("Counts base statistics in FASTA file(s)")
        .after_help(
            r###"
This command calculates the base statistics (A, C, G, T, N) for each sequence in one or more FASTA files.
It outputs a TSV table with the following columns:
* seq: Sequence name
* len: Sequence length
* A, C, G, T, N: Count of each base
* ignored: Count of other characters (e.g., IUPAC codes, gaps)

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

Examples:
1. Count base statistics for a single FASTA file:
   pgr fa count input.fa

2. Count base statistics for multiple FASTA files:
   pgr fa count input1.fa input2.fa
"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input FASTA file(s) to process"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    //----------------------------
    // Init
    //----------------------------
    let mut total_len = 0usize;
    let mut total_base_cnt = [0usize; 5]; // A, C, G, T, N

    // Write the header
    writer.write_fmt(format_args!("#seq\tlen\tA\tC\tG\tT\tN\n"))?;

    //----------------------------
    // Process
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut fa_in = pgr::libs::fmt::fa::reader(infile)?;

        for result in fa_in.records() {
            let record = result?;
            let name = String::from_utf8(record.name().into())?;
            let seq = record.sequence();

            let (len, base_cnt) = pgr::libs::fasta::stat::count_bases(seq.get(..).unwrap());

            writer.write_fmt(format_args!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                name,
                len,
                base_cnt[pgr::libs::nt::Nt::A as usize],
                base_cnt[pgr::libs::nt::Nt::C as usize],
                base_cnt[pgr::libs::nt::Nt::G as usize],
                base_cnt[pgr::libs::nt::Nt::T as usize],
                base_cnt[pgr::libs::nt::Nt::N as usize],
            ))?;

            // Update total statistics
            total_len += len;
            for &nt in &[
                pgr::libs::nt::Nt::A,
                pgr::libs::nt::Nt::C,
                pgr::libs::nt::Nt::G,
                pgr::libs::nt::Nt::T,
                pgr::libs::nt::Nt::N,
            ] {
                total_base_cnt[nt as usize] += base_cnt[nt as usize];
            }
        }
    }

    //----------------------------
    // Output total
    //----------------------------
    writer.write_fmt(format_args!(
        "total\t{}\t{}\t{}\t{}\t{}\t{}\n",
        total_len,
        total_base_cnt[pgr::libs::nt::Nt::A as usize],
        total_base_cnt[pgr::libs::nt::Nt::C as usize],
        total_base_cnt[pgr::libs::nt::Nt::G as usize],
        total_base_cnt[pgr::libs::nt::Nt::T as usize],
        total_base_cnt[pgr::libs::nt::Nt::N as usize],
    ))?;

    Ok(())
}
