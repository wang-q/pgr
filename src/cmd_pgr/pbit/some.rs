use anyhow::Context;
use clap::{ArgMatches, Command};
use pgr::libs::pbit::decompressor::Decompressor;
use std::io::Write;

/// Build the clap subcommand for some.
pub fn make_subcommand() -> Command {
    Command::new("some")
        .about("Extracts sample sequences based on a list of contig names")
        .after_help(
            r###"
This command extracts full contig sequences from all samples in a pbit
archive, filtered by a list of contig names.

Notes:
* Case-sensitive name matching
* One contig name per line in the list file
* Empty lines and lines starting with '#' are ignored
* Output format is FASTA (one entry per sample-contig pair)
* pbit files are binary and require random access (seeking)
* Does not support stdin or gzipped inputs

Examples:
1. Extract contigs listed in list.txt:
   pgr pbit some input.pbit list.txt -o output.fa

2. Extract contigs NOT in list.txt:
   pgr pbit some input.pbit list.txt -i -o output.fa
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input pbit file to process",
        ))
        .arg(crate::cmd_pgr::args::fa_name_list_arg(true))
        .arg(crate::cmd_pgr::args::invert_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the some command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let is_invert = args.get_flag("invert");
    let infile = args.get_one::<String>("infile").unwrap();
    let list_file = args.get_one::<String>("name_list").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    // Load contig name list.
    let set_list = pgr::libs::io::read_names::<std::collections::HashSet<String>>(list_file)?;

    let mut dec = Decompressor::open(infile)
        .with_context(|| format!("Failed to open pbit file {}", infile))?;
    let mut writer = pgr::libs::io::writer(outfile)
        .with_context(|| format!("Failed to open writer for {}", outfile))?;

    // Gather all contig names across all samples (deduplicated, ordered).
    // Collect into owned Strings to release the immutable borrow on dec
    // before calling dec.get_contig (mutable).
    let all_contigs: Vec<String> = dec
        .list_contigs(None)
        .into_iter()
        .map(String::from)
        .collect();

    for contig in &all_contigs {
        let keep = set_list.contains(contig.as_str()) != is_invert;
        if !keep {
            continue;
        }
        // Extract full contig from all samples (no slice, positive strand).
        dec.get_contig(contig, None, None, "+", &mut writer)?;
    }

    writer.flush()?;
    Ok(())
}
