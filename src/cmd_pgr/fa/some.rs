use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("some")
        .about("Extracts FASTA records based on a list of names")
        .after_help(
            r###"
This command extracts FASTA records from an input file based on a list of sequence names.

Notes:
* Case-sensitive name matching
* One sequence name per line in the list file
* Empty lines and lines starting with '#' are ignored
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

Examples:
1. Extract sequences listed in list.txt:
   pgr fa some input.fa list.txt

2. Extract sequences NOT in list.txt:
   pgr fa some input.fa list.txt -i

3. Process gzipped files:
   pgr fa some input.fa.gz list.txt -o output.fa.gz

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input FASTA file to process"),
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
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let is_invert = args.get_flag("invert");

    let mut fa_in = pgr::libs::fmt::fa::reader(args.get_one::<String>("infile").unwrap())?;

    let mut fa_out = pgr::libs::fmt::fa::writer(crate::cmd_pgr::args::get_outfile(args))?;

    //----------------------------
    // Load list
    //----------------------------
    let set_list = pgr::libs::io::read_names_as_set(args.get_one::<String>("list.txt").unwrap())?;

    //----------------------------
    // Process
    //----------------------------
    for result in fa_in.records() {
        let record = result?;
        let name = String::from_utf8(record.name().into())?;

        if set_list.contains(&name) != is_invert {
            fa_out.write_record(&record)?;
        }
    }

    Ok(())
}
