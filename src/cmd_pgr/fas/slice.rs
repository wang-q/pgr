use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for slice.
pub fn make_subcommand() -> Command {
    Command::new("slice")
        .about("Extracts alignment slices")
        .after_help(
            r###"
Extracts alignment slices from block FA files using a runlist JSON.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* The JSON file (--runlist) keys are chromosome/sequence names, and values are runlists (e.g., "1-100,200-300")
* If `--name` is not specified, the first species of the first non-empty block is used as the reference

Examples:
1. Extract slices defined in a JSON file:
   pgr fas slice tests/fas/slice.fas --runlist tests/fas/slice.json

2. Extract slices and name the output based on a specific species:
   pgr fas slice tests/fas/slice.fas --runlist tests/fas/slice.json --name S288c

3. Output results to a file:
   pgr fas slice tests/fas/slice.fas --runlist tests/fas/slice.json -o output.fas

"###,
        )
        .arg(crate::cmd_pgr::args::runlist_arg())
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(
            crate::cmd_pgr::args::fas_name_arg("Reference species name. Default is the first species"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the slice command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let runlist = args.get_one::<String>("runlist").unwrap();

    // Protect both the block FA inputs and the runlist JSON: the writer is
    // opened (truncating) before any of them is read.
    let mut inputs: Vec<&str> = args
        .get_many::<String>("infiles")
        .unwrap()
        .map(|s| s.as_str())
        .collect();
    inputs.push(runlist.as_str());
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, inputs)?;

    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    let set = pgr::libs::io::read_runlist(runlist)?;

    let mut name = args
        .get_one::<String>("name")
        .map(|s| s.to_string())
        .unwrap_or_default();

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader =
            pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;

        for block_result in pgr::libs::fmt::fas::iter_fas_blocks(&mut reader) {
            let block = block_result?;
            if block.entries.is_empty() {
                continue;
            }
            // the first name of the first non-empty block becomes the default reference
            if name.is_empty() {
                name = block.names[0].clone();
            }

            pgr::libs::alignment::slice_block(&block, &name, &set, &mut writer)?;
        }
    }

    writer.flush()?;
    Ok(())
}
