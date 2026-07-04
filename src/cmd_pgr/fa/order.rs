use clap::{ArgMatches, Command};
use std::collections::BTreeMap;

/// Build the clap subcommand for order.
pub fn make_subcommand() -> Command {
    Command::new("order")
        .about("Extracts FASTA records in the order specified by a list")
        .after_help(
            r###"
This command extracts FASTA records from an input file in the order specified by a list of sequence names.

Notes:
* Case-sensitive name matching
* One sequence name per line in the list file
* Empty lines and lines starting with '#' are ignored
* All sequences are loaded into memory
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* Missing sequences in the input file are silently skipped

Examples:
1. Extract sequences in order specified by list.txt:
   pgr fa order input.fa list.txt

2. Process gzipped files:
   pgr fa order input.fa.gz list.txt -o output.fa.gz

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input FASTA file to process",
        ))
        .arg(crate::cmd_pgr::args::fa_name_list_arg(true))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the order command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut fa_in = pgr::libs::fmt::fa::reader(args.get_one::<String>("infile").unwrap())?;

    let mut fa_out = pgr::libs::fmt::fa::writer(crate::cmd_pgr::args::get_outfile(args))?;

    let list: indexmap::IndexSet<_> =
        pgr::libs::io::read_names::<Vec<String>>(args.get_one::<String>("name_list").unwrap())?
            .into_iter()
            .collect();

    // Load records into a BTreeMap for efficient lookup
    let mut record_of = BTreeMap::new();

    for result in fa_in.records() {
        let record = result?;
        let name = String::from_utf8(record.name().into())?;

        if list.contains(&name) {
            record_of.insert(name, record);
        }
    }

    for name in list.iter() {
        if let Some(record) = record_of.get(name) {
            fa_out.write_record(record)?;
        }
    }

    Ok(())
}
