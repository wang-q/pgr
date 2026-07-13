use anyhow::Context;
use clap::{ArgMatches, Command};
use indexmap::IndexMap;
use std::io::Write;

/// Build the clap subcommand for name.
pub fn make_subcommand() -> Command {
    Command::new("name")
        .about("Outputs all species names from block FA files")
        .after_help(
            r###"
Extracts and outputs all species names from block FA files.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* By default, the subcommand outputs a list of unique species names
* Use `--count` to also output the number of occurrences of each species name

Examples:
1. Output all species names:
   pgr fas name tests/fas/example.fas

2. Output species names with occurrence counts:
   pgr fas name tests/fas/example.fas --count

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(crate::cmd_pgr::args::count_arg(
            "Output species names with occurrence counts",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the name command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    let is_count = args.get_flag("count");

    // Operating
    let mut counts: IndexMap<String, i32> = IndexMap::new();

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader =
            pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;

        for block_result in pgr::libs::fmt::fas::iter_fas_blocks(&mut reader) {
            let block = block_result?;
            for entry in &block.entries {
                let name = entry.range().name().to_string();
                *counts.entry(name).or_insert(0) += 1;
            }
        }
    }

    for (name, count) in &counts {
        if is_count {
            writer.write_all(format!("{}\t{}\n", name, count).as_ref())?;
        } else {
            writer.write_all(format!("{}\n", name).as_ref())?;
        }
    }

    writer.flush()?;
    Ok(())
}
