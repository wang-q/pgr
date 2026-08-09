use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for six-frame.
pub fn make_subcommand() -> Command {
    Command::new("six-frame")
        .about("Translates DNA sequences in six frames")
        .after_help(
            r###"
This command translates DNA sequences in six frames and identifies ORFs.

Output format:
>sequence_name(strand):start-end|frame=N
MXXXXXX*

Translation frames:
* frame is the 0-based reading-frame offset within the strand (0, 1, or 2)
* Forward strand: strand is '+'; the offset is applied to the sequence as-is
* Reverse strand: strand is '-'; the offset is applied after reverse-complementing

Notes:
* Filters: --min-len N (min length), --start-met (starts with M), --end (ends with *)
* Coordinates are 1-based
* Non-standard bases are translated as X
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* Stop codons are included in the output

Examples:
1. Basic translation:
   pgr fa six-frame input.fa -o orfs.fa

2. Filter long ORFs:
   pgr fa six-frame input.fa --min-len 100 -o orfs.fa

3. Complete proteins only:
   pgr fa six-frame input.fa --start-met --end -o orfs.fa

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input FASTA file to process",
        ))
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

/// Execute the six-frame command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let mut reader = pgr::libs::fmt::seq::SeqReader::new(infile)
        .with_context(|| format!("Failed to open reader for {}", infile))?;
    let mut rec = pgr::libs::fmt::seq::SeqRecord::new();

    let opt_len = *args.get_one::<usize>("min_len").unwrap();
    let is_start = args.get_flag("start_met");
    let is_end = args.get_flag("end");

    let outfile = crate::cmd_pgr::args::get_outfile(args);
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, [infile.as_str()])?;
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    while reader.read_record(&mut rec)? {
        let name = String::from_utf8(rec.name().to_vec())?;
        let seq = rec.sequence();

        // Perform six-frame translation
        let translations = pgr::libs::translate::six_frame_translation(seq);

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

    writer.flush()?;
    Ok(())
}
