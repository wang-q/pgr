use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("one")
        .about("Extracts one FASTA record by name")
        .after_help(
            r###"
This command extracts a single record from a FASTA file by its sequence name (ID).

Notes:
* Scans the file sequentially to find the matching record
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

Examples:
1. Extract a record by name:
   pgr fa one input.fa chr1
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input FASTA file to process",
        ))
        .arg(
            Arg::new("name")
                .required(true)
                .index(2)
                .help("Name of the sequence to extract"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let mut fa_in = pgr::libs::fmt::fa::reader(args.get_one::<String>("infile").unwrap())?;

    let mut fa_out = pgr::libs::fmt::fa::writer(crate::cmd_pgr::args::get_outfile(args))?;

    let name = args.get_one::<String>("name").unwrap();

    //----------------------------
    // Process
    //----------------------------
    for result in fa_in.records() {
        let record = result?;
        let this_name = String::from_utf8(record.name().into())?;

        if this_name == *name {
            fa_out.write_record(&record)?;
            break;
        }
    }

    Ok(())
}
