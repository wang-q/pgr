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
* Reference genome can be plain text or bgzipped
* Output format: `range<TAB>status` where status is OK or FAILED

Examples:"###,
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
    let opt_genome = args.get_one::<String>("genome").unwrap();
    let opt_name: &str = args
        .get_one::<String>("name")
        .map(|s| s.as_str())
        .unwrap_or("");

    // The writer is opened (truncating) before the block FA inputs and the
    // reference genome are read, so reject an `-o` that would overwrite any of
    // them. The genome's `.loc` sidecar is also protected: truncating it would
    // make `open_indexed` treat it as fresh and serve an empty index.
    let mut inputs: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .map(|s| s.to_string())
        .collect();
    inputs.push(opt_genome.to_string());
    inputs.push(format!("{}.loc", opt_genome));
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, inputs.iter().map(|s| s.as_str()))?;

    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

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
