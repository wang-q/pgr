use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("count")
        .about("Counts base statistics in FASTA file(s)")
        .after_help(
            r###"
This command calculates the base statistics (A, C, G, T, N) for each sequence in one or more FASTA files.
It outputs a TSV table with the following columns:
1. seq: Sequence name
2. len: Sequence length
3. A, C, G, T, N: Count of each base
4. ignored: Count of other characters (e.g., IUPAC codes, gaps)

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
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

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
        let reader = intspan::reader(infile);
        let mut fa_in = noodles_fasta::io::Reader::new(reader);

        for result in fa_in.records() {
            let record = result?;
            let name = String::from_utf8(record.name().into())?;
            let seq = record.sequence();

            let (len, base_cnt) = count_bases(seq.get(..).unwrap());

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

// Count bases in a sequence
fn count_bases(seq: &[u8]) -> (usize, [usize; 5]) {
    let mut len = 0usize;
    let mut base_cnt = [0usize; 5]; // A, C, G, T, N

    for &el in seq {
        let nt = pgr::libs::nt::to_nt(el);
        if !matches!(nt, pgr::libs::nt::Nt::Invalid) {
            len += 1;
            base_cnt[nt as usize] += 1;
        }
    }

    (len, base_cnt)
}
