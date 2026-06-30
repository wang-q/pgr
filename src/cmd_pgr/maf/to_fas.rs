use clap::*;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("to-fas")
        .about("Convert MAF files to block FA format")
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
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input MAF file(s) to process"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    //----------------------------
    // Ops
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::reader(infile)?;

        while let Ok(block) = pgr::libs::fmt::maf::next_maf_block(&mut reader) {
            for comp in block.components {
                let range = comp.to_range();

                //----------------------------
                // Output
                //----------------------------
                writer.write_all(format!(">{}\n{}\n", range, comp.text).as_ref())?;
            }

            // end of a block
            writer.write_all("\n".as_ref())?;
        }
    }

    Ok(())
}
