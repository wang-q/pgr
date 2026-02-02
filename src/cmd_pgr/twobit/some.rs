use clap::*;
use pgr::libs::twobit::TwoBitFile;
use std::collections::HashSet;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("some")
        .about("Extract some records from a 2bit file based on a list of names")
        .after_help(
            r###"
This command extracts sequences from a 2bit file based on a list of sequence names
and writes them to a FASTA file.

Notes:
* Case-sensitive name matching
* One sequence name per line in the list file
* Empty lines and lines starting with '#' are ignored
* Input must be a 2bit file (stdin not supported due to seeking requirements)
* Output format is FASTA

Examples:
1. Extract sequences listed in list.txt:
   pgr 2bit some input.2bit list.txt -o output.fa

2. Extract sequences NOT in list.txt:
   pgr 2bit some input.2bit list.txt -i -o output.fa

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input 2bit file to process"),
        )
        .arg(
            Arg::new("list.txt")
                .required(true)
                .index(2)
                .help("File containing one sequence name per line"),
        )
        .arg(
            Arg::new("invert")
                .long("invert")
                .short('i')
                .action(ArgAction::SetTrue)
                .help("Invert selection: output sequences NOT in the list"),
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
    let is_invert = args.get_flag("invert");
    let infile = args.get_one::<String>("infile").unwrap();
    let list_file = args.get_one::<String>("list.txt").unwrap();
    let outfile = args.get_one::<String>("outfile").unwrap();

    //----------------------------
    // Load list
    //----------------------------
    let set_list: HashSet<String> = intspan::read_first_column(list_file)
        .into_iter()
        .collect();

    //----------------------------
    // Process
    //----------------------------
    let mut tb = TwoBitFile::open(infile)?;
    let names = tb.get_sequence_names();

    let mut writer = intspan::writer(outfile);

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
