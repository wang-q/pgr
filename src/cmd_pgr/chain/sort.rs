use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::chain::{read_chains, ChainReader};
use std::io::BufRead;
use std::io::Write;
/// Build the clap subcommand for sort.
pub fn make_subcommand() -> Command {
    Command::new("sort")
        .about("Sorts chains by score")
        .after_help(
            r###"
Sorts chains by score in descending order. By default, chain IDs are renumbered
starting from 1 after sorting; use `--save-id` to preserve the original IDs.

Notes:
* Accepts multiple input files; they are concatenated then sorted together
* Use `--input-list` to read input file paths from a list (one per line)
* Output is written to stdout if `--outfile` is omitted

Examples:
1. Sort a single chain file:
   pgr chain sort in.chain -o sorted.chain

2. Preserve original chain IDs:
   pgr chain sort in.chain --save-id -o sorted.chain

3. Concatenate and sort from a file list:
   pgr chain sort --input-list files.txt -o sorted.chain

"###,
        )
        .arg(
            Arg::new("infiles")
                .required_unless_present("input_list")
                .num_args(1..)
                .action(ArgAction::Append)
                .help("Input chain file(s)"),
        )
        .arg(
            Arg::new("input_list")
                .long("input-list")
                .help("File containing a list of input chain files (one per line)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("save_id")
                .long("save-id")
                .action(ArgAction::SetTrue)
                .help("Keep existing chain IDs (default: renumber starting from 1)"),
        )
}
/// Execute the sort command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut files: Vec<String> = args
        .get_many::<String>("infiles")
        .map(|v| v.cloned().collect())
        .unwrap_or_default();

    if let Some(list_path) = args.get_one::<String>("input_list") {
        let reader = pgr::reader(list_path)
            .with_context(|| format!("Failed to open reader for {}", list_path))?;
        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if !trimmed.is_empty() {
                files.push(trimmed.to_string());
            }
        }
    }

    let save_id = args.get_flag("save_id");

    let mut all_chains = Vec::new();
    let mut header_comments: Vec<String> = Vec::new();

    // Read all chains. UCSC chainMergeSort propagates `##` metadata from the
    // first file only (lineFileSetMetaDataOutput on the first lineFile).
    for (i, file_path) in files.iter().enumerate() {
        let reader = pgr::reader(file_path)
            .with_context(|| format!("Failed to open reader for {}", file_path))?;
        if i == 0 {
            let mut chain_reader = ChainReader::new(reader);
            all_chains.extend(chain_reader.by_ref().collect::<Result<Vec<_>, _>>()?);
            header_comments = std::mem::take(&mut chain_reader.header_comments);
        } else {
            let chains = read_chains(reader)?;
            all_chains.extend(chains);
        }
    }

    // Sort by score descending, renumber unless --save-id
    pgr::libs::chain::sort_chains(&mut all_chains, !save_id);

    // Write output
    let out_path = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(out_path).with_context(|| format!("Failed to open writer for {}", out_path))?;
    for comment in &header_comments {
        writeln!(writer, "{}", comment.trim_end())?;
    }
    for chain in all_chains {
        chain.write(&mut writer)?;
    }

    writer.flush()?;
    Ok(())
}
