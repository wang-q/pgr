use anyhow::Context;
use clap::{ArgMatches, Command};
use std::collections::BTreeMap;
use std::io::Write;

/// Build the clap subcommand for join.
pub fn make_subcommand() -> Command {
    Command::new("join")
        .about("Joins multiple block fasta files by a common target")
        .after_help(
            r###"
Joins multiple block fasta files by a common target sequence.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

Examples:
1. Join multiple block FA files:
   pgr fas join tests/fas/part1.fas tests/fas/part2.fas

2. Join files based on a specific species:
   pgr fas join tests/fas/part1.fas tests/fas/part2.fas --name S288c

3. Output results to a file:
   pgr fas join tests/fas/part1.fas tests/fas/part2.fas -o output.fas

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(crate::cmd_pgr::args::fas_name_arg(
            "Target species name. Default is the first species",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the join command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    let mut name = if args.contains_id("name") {
        args.get_one::<String>("name").unwrap().to_string()
    } else {
        "".to_string()
    };
    let mut block_of: BTreeMap<String, Vec<pgr::libs::fmt::fas::FasEntry>> = BTreeMap::new();

    // Operating
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader =
            pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;

        for block_result in pgr::libs::fmt::fas::iter_fas_blocks(&mut reader) {
            let block = block_result?;
            if name.is_empty() {
                name = block
                    .names
                    .first()
                    .ok_or_else(|| anyhow::anyhow!("empty block names"))?
                    .to_string();
            }

            pgr::libs::fmt::fas::join_block_entries(&block, &name, &mut block_of)?;
        }
    }

    for v in block_of.values() {
        for e in v {
            writer.write_all(e.to_string().as_ref())?;
        }
        // end of a block
        writer.write_all("\n".as_ref())?;
    }

    writer.flush()?;
    Ok(())
}
