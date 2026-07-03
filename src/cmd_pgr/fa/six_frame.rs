use clap::*;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("six-frame")
        .about("Translates DNA sequences in six frames")
        .after_help(
            r###"
This command translates DNA sequences in six frames and identifies ORFs.

Translation frames:
* Forward strand: +1, +2, +3 (starting at positions 0, 1, 2)
* Reverse strand: -1, -2, -3 (complement sequence, then start at 0, 1, 2)

Output format:
>sequence_name(strand):start-end|frame=N
MXXXXXX*

Notes:
* Filters: --len N (min length), --start (starts with M), --end (ends with *)
* Coordinates are 1-based
* Non-standard bases are translated as X
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* Stop codons are included in the output

Examples:
1. Basic translation:
   pgr fa six-frame input.fa -o orfs.fa

2. Filter long ORFs:
   pgr fa six-frame input.fa --len 100 -o orfs.fa

3. Complete proteins only:
   pgr fa six-frame input.fa --start --end -o orfs.fa

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Input FASTA file to process"),
        )
        .arg(crate::cmd_pgr::args::min_len_arg_with_default(
            "0",
            "Minimum length of the amino acid sequence to consider",
        ))
        .arg(
            Arg::new("start_met")
                .long("start-met")
                .action(ArgAction::SetTrue)
                .help("Only consider ORFs that start with Methionine (M)"),
        )
        .arg(
            Arg::new("end")
                .long("end")
                .action(ArgAction::SetTrue)
                .help("Only consider ORFs that end with a stop codon (*)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let mut fa_in = pgr::libs::fmt::fa::reader(args.get_one::<String>("infile").unwrap())?;

    let opt_len = *args.get_one::<usize>("min_len").unwrap();
    let is_start = args.get_flag("start_met");
    let is_end = args.get_flag("end");

    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    //----------------------------
    // Ops
    //----------------------------
    for result in fa_in.records() {
        // obtain record or fail with error
        let record = result?;

        let name = String::from_utf8(record.name().into())?;
        let seq = record.sequence();

        // Perform six-frame translation
        let translations = pgr::libs::translate::six_frame_translation(&seq[..]);

        for (protein, frame, is_reverse) in translations {
            let orfs = pgr::libs::translate::find_orfs(&protein);
            let filtered = pgr::libs::translate::filter_and_convert_orfs(
                &orfs,
                seq.len(),
                frame,
                is_reverse,
                opt_len,
                is_start,
                is_end,
            );

            for (orf_start, orf_end, orf_seq) in filtered {
                let header = format!(
                    "{}({}):{}-{}|frame={}",
                    name,
                    if is_reverse { "-" } else { "+" },
                    orf_start,
                    orf_end,
                    frame,
                );
                writer.write_fmt(format_args!(">{}\n{}\n", header, orf_seq))?;
            }
        }
    }

    Ok(())
}
