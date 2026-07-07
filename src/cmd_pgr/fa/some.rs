use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for some.
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
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input FASTA file to process",
        ))
        .arg(crate::cmd_pgr::args::fa_name_list_arg(true))
        .arg(crate::cmd_pgr::args::invert_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the some command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let is_invert = args.get_flag("invert");

    let infile = args.get_one::<String>("infile").unwrap();
    let mut fa_in = pgr::libs::fmt::fa::reader(infile)
        .with_context(|| format!("Failed to open reader for {}", infile))?;

    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut fa_out = pgr::libs::fmt::fa::writer(outfile)
        .with_context(|| format!("Failed to open writer for {}", outfile))?;

    // Load list
    let set_list = pgr::libs::io::read_names::<std::collections::HashSet<String>>(
        args.get_one::<String>("name_list").unwrap(),
    )?;

    for result in fa_in.records() {
        let record = result?;
        let name = String::from_utf8(record.name().into())?;

        if set_list.contains(&name) != is_invert {
            fa_out.write_record(&record)?;
        }
    }

    fa_out.get_mut().flush()?;

    Ok(())
}
