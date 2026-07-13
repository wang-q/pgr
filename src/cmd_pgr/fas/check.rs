use anyhow::Context;
use clap::{ArgMatches, Command};
use pgr::libs::loc;
use std::io::Write;

/// Build the clap subcommand for check.
pub fn make_subcommand() -> Command {
    Command::new("check")
        .about("Checks genome locations in block FA headers")
        .after_help(
            r###"
Checks genome locations in block FA headers against a reference genome FA file.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

"###,
        )
        .arg(crate::cmd_pgr::args::genome_arg())
        .arg(crate::cmd_pgr::args::infiles_arg_with_help(
            "Input block FA file(s) to check",
        ))
        .arg(crate::cmd_pgr::args::fas_name_arg(
            "Check sequences for a specific species",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the check command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    let opt_genome = args.get_one::<String>("genome").unwrap();
    let opt_name: &str = args
        .get_one::<String>("name")
        .map(|s| s.as_str())
        .unwrap_or("");

    let (mut genome_reader, loc_of) = loc::open_indexed(opt_genome, false)?;

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader =
            pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;

        for block_result in pgr::libs::fmt::fas::iter_fas_blocks(&mut reader) {
            let block = block_result?;

            for (entry, name) in block.entries.iter().zip(&block.names) {
                if !opt_name.is_empty() && name != opt_name {
                    continue;
                }
                let status = pgr::libs::fmt::fas::check_entry_against_ref(
                    entry,
                    &mut genome_reader,
                    &loc_of,
                )?;
                writer.write_all(format!("{}\t{}\n", entry.range(), status).as_ref())?;
            }
        }
    }

    writer.flush()?;
    Ok(())
}
