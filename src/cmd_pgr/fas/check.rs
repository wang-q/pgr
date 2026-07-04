use clap::{ArgMatches, Command};
use pgr::libs::loc;
use std::io::Write;

/// Build the clap subcommand for check.
pub fn make_subcommand() -> Command {
    Command::new("check")
        .about("Checks genome locations in block FA headers")
        .after_help(
            r###"
Checks genome locations in block FA headers against a chrom.sizes file.

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
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;
    let opt_genome = args.get_one::<String>("genome").unwrap();
    let opt_name = &args
        .get_one::<String>("name")
        .map(|s| s.as_str())
        .unwrap_or("")
        .to_string();

    let (mut genome_reader, loc_of) = loc::open_indexed(opt_genome, false)?;

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::reader(infile)?;

        for block_result in pgr::libs::fmt::fas::iter_fas_blocks(&mut reader) {
            let block = block_result?;
            let block_names = &block.names;

            // Check if a specific species is requested
            if !opt_name.is_empty() && block_names.contains(opt_name) {
                for entry in &block.entries {
                    let entry_name = entry.range().name();
                    if entry_name == opt_name {
                        let status = pgr::libs::fmt::fas::check_entry_against_ref(
                            entry,
                            &mut genome_reader,
                            &loc_of,
                        )?;
                        writer.write_all(format!("{}\t{}\n", entry.range(), status).as_ref())?;
                    }
                }
            } else if opt_name.is_empty() {
                // Check all sequences in the block
                for entry in &block.entries {
                    let status = pgr::libs::fmt::fas::check_entry_against_ref(
                        entry,
                        &mut genome_reader,
                        &loc_of,
                    )?;
                    writer.write_all(format!("{}\t{}\n", entry.range(), status).as_ref())?;
                }
            }
        }
    }

    Ok(())
}
