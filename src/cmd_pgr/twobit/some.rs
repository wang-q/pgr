use clap::*;
use pgr::libs::fmt::twobit::TwoBitFile;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("some")
        .about("Extracts 2bit records based on a list of names")
        .after_help(
            r###"
This command extracts sequences from a 2bit file based on a list of sequence names
and writes them to a FASTA file.

Notes:
* Case-sensitive name matching
* One sequence name per line in the list file
* Empty lines and lines starting with '#' are ignored
* Output format is FASTA
* 2bit files are binary and require random access (seeking)
* Does not support stdin or gzipped inputs

Examples:
1. Extract sequences listed in list.txt:
   pgr 2bit some input.2bit list.txt -o output.fa

2. Extract sequences NOT in list.txt:
   pgr 2bit some input.2bit list.txt -i -o output.fa

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input 2bit file to process",
        ))
        .arg(
            Arg::new("list")
                .required(true)
                .index(2)
                .help("File containing one sequence name per line"),
        )
        .arg(crate::cmd_pgr::args::invert_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let is_invert = args.get_flag("invert");
    let infile = args.get_one::<String>("infile").unwrap();
    let list_file = args.get_one::<String>("list").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    //----------------------------
    // Load list
    //----------------------------
    let set_list = pgr::libs::io::read_names::<std::collections::HashSet<String>>(list_file)?;

    //----------------------------
    // Process
    //----------------------------
    let mut tb = TwoBitFile::open(infile)?;
    let names = tb.get_sequence_names();

    let mut writer = pgr::writer(outfile)?;

    for name in names {
        if set_list.contains(&name) != is_invert {
            // Read sequence with masking (no_mask = false)
            let seq = tb.read_sequence(&name, None, None, false)?;

            // Write FASTA
            // Matches pgr fa some behavior (single line sequence)
            write!(writer, ">{}\n{}\n", name, seq)?;
        }
    }

    Ok(())
}
