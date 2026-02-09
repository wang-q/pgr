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
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap());

    //----------------------------
    // Ops
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::reader(infile);

        while let Ok(block) = pgr::libs::fas::next_maf_block(&mut reader) {
            // Can't use reference as entry.alignment does not Copy
            for entry in block.entries {
                let range = entry.to_range();
                let seq = String::from_utf8(entry.alignment).unwrap();

                //----------------------------
                // Output
                //----------------------------
                writer.write_all(format!(">{}\n{}\n", range, seq).as_ref())?;
            }

            // end of a block
            writer.write_all("\n".as_ref())?;
        }
    }

    Ok(())
}
