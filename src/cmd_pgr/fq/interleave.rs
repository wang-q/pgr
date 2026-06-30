use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("interleave")
        .visible_alias("il")
        .about("Interleave paired-end sequences")
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
   pgr fq interleave R1.fa --prefix sample --start 1

3. Convert to FQ:
   pgr fq interleave R1.fa R2.fa --fq -o out.fq

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..=2)
                .index(1)
                .help("Set the input files to use"),
        )
        .arg(
            Arg::new("fq")
                .long("fq")
                .action(ArgAction::SetTrue)
                .help("Write FQ"),
        )
        .arg(
            Arg::new("prefix")
                .long("prefix")
                .num_args(1)
                .default_value("read")
                .help("Prefix of record names"),
        )
        .arg(
            Arg::new("start")
                .long("start")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("0")
                .help("Starting index"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap())?;

    let is_out_fq = args.get_flag("fq");
    let opt_prefix = args.get_one::<String>("prefix").unwrap();
    let mut opt_start = *args.get_one::<usize>("start").unwrap();

    let infiles: Vec<&str> = args
        .get_many::<String>("infiles")
        .unwrap()
        .map(|s| s.as_str())
        .collect();
    let is_in_fq = pgr::libs::io::is_fq(infiles[0])?;

    if infiles.len() == 1 {
        // single file: produce dummy R2 for each R1
        if is_in_fq {
            let reader = pgr::reader(infiles[0])?;
            let mut seq_in = noodles_fastq::io::Reader::new(reader);
            for result in seq_in.records() {
                let record = result?;
                // Preserve original: dummy R2 seq is "\n" for FA output, "N" for FQ output
                let r2_seq: &[u8] = if is_out_fq { b"N" } else { b"\n" };
                pgr::libs::fmt::fq::write_pair(
                    &mut writer,
                    opt_prefix,
                    opt_start,
                    record.sequence(),
                    Some(record.quality_scores()),
                    r2_seq,
                    Some(b"!"),
                    is_out_fq,
                )?;
                opt_start += 1;
            }
        } else {
            let mut seq_in = pgr::libs::fmt::fa::reader(infiles[0])?;
            for result in seq_in.records() {
                let record = result?;
                pgr::libs::fmt::fq::write_pair(
                    &mut writer,
                    opt_prefix,
                    opt_start,
                    &record.sequence()[..],
                    None,
                    b"N",
                    None,
                    is_out_fq,
                )?;
                opt_start += 1;
            }
        }
    } else {
        // two files: zip R1 and R2 records
        if is_in_fq {
            let reader1 = pgr::reader(infiles[0])?;
            let mut seq1_in = noodles_fastq::io::Reader::new(reader1);
            let reader2 = pgr::reader(infiles[1])?;
            let mut seq2_in = noodles_fastq::io::Reader::new(reader2);
            for (r1, r2) in std::iter::zip(seq1_in.records(), seq2_in.records()) {
                let record1 = r1?;
                let record2 = r2?;
                pgr::libs::fmt::fq::write_pair(
                    &mut writer,
                    opt_prefix,
                    opt_start,
                    record1.sequence(),
                    Some(record1.quality_scores()),
                    record2.sequence(),
                    Some(record2.quality_scores()),
                    is_out_fq,
                )?;
                opt_start += 1;
            }
        } else {
            let mut seq1_in = pgr::libs::fmt::fa::reader(infiles[0])?;
            let mut seq2_in = pgr::libs::fmt::fa::reader(infiles[1])?;
            for (r1, r2) in std::iter::zip(seq1_in.records(), seq2_in.records()) {
                let record1 = r1?;
                let record2 = r2?;
                pgr::libs::fmt::fq::write_pair(
                    &mut writer,
                    opt_prefix,
                    opt_start,
                    &record1.sequence()[..],
                    None,
                    &record2.sequence()[..],
                    None,
                    is_out_fq,
                )?;
                opt_start += 1;
            }
        }
    }

    Ok(())
}
