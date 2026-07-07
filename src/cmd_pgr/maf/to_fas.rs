use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for to-fas.
pub fn make_subcommand() -> Command {
    Command::new("to-fas")
        .about("Converts MAF files to block FA format")
        .after_help(
            r###"
This subcommand converts MAF (Multiple Alignment Format) files into block FA format.

Input files can be gzipped. If the input file is 'stdin', data is read from standard input.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* MAF files typically contain multiple sequence alignments, and this tool extracts each alignment into block FASTA format
* The output preserves the alignment structure, with each block separated by a newline

Examples:
1. Convert a MAF file to block FASTA format:
   pgr maf to-fas tests/maf/example.maf

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("MAF"))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the to-fas command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader =
            pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;

        loop {
            let block = match pgr::libs::fmt::maf::next_maf_block(&mut reader) {
                Ok(b) => b,
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(e.into()),
            };
            for comp in block.components {
                let range = comp.to_range();

                writer.write_all(format!(">{}\n{}\n", range, comp.text).as_ref())?;
            }

            // end of a block
            writer.write_all("\n".as_ref())?;
        }
    }

    writer.flush()?;
    Ok(())
}
