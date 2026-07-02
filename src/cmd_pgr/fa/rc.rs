use clap::*;
use std::collections::HashSet;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("rc")
        .about("Reverse complements sequences in FASTA file(s)")
        .after_help(
            r###"
This command reverse complements DNA sequences in FASTA files.

Notes:
* Process all sequences or only selected ones
* Optionally prefix names with 'RC_'
* Handles IUPAC ambiguous codes correctly
* Preserves case (upper/lower) of bases
* Case-sensitive name matching when using list
* Empty lines and lines starting with '#' are ignored in list
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* Non-IUPAC characters are preserved as-is

Examples:
1. Reverse complement all sequences:
   pgr fa rc input.fa -o output.fa

2. Only process listed sequences:
   pgr fa rc input.fa list.txt -o output.fa

3. Keep original names (no 'RC_' prefix):
   pgr fa rc input.fa -c -o output.fa

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input FASTA file to process"),
        )
        .arg(
            Arg::new("list")
                .required(false)
                .index(2)
                .help("File containing one sequence name per line (optional)"),
        )
        .arg(
            Arg::new("consistent")
                .long("consistent")
                .short('c')
                .action(ArgAction::SetTrue)
                .help("Keep the name consistent (don't prepend 'RC_')"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let mut fa_in = pgr::libs::fmt::fa::reader(args.get_one::<String>("infile").unwrap())?;

    let is_consistent = args.get_flag("consistent");

    let mut fa_out = pgr::libs::fmt::fa::writer(crate::cmd_pgr::args::get_outfile(args))?;

    let set_list: HashSet<String> = if args.contains_id("list") {
        pgr::libs::io::read_names::<std::collections::HashSet<String>>(
            args.get_one::<String>("list").unwrap(),
        )?
    } else {
        HashSet::new()
    };

    //----------------------------
    // Process
    //----------------------------
    for result in fa_in.records() {
        let record = result?;
        let name = String::from_utf8(record.name().into())?;

        if args.contains_id("list") && !set_list.contains(&name) {
            fa_out.write_record(&record)?;
            continue;
        }

        let new_name = if is_consistent {
            name
        } else {
            format!("RC_{}", name)
        };

        let seq_rc: Vec<u8> = record
            .sequence()
            .complement()
            .rev()
            .collect::<Result<_, _>>()?;
        let record_rc = pgr::libs::fmt::fa::new_record(&new_name, &seq_rc);
        fa_out.write_record(&record_rc)?;
    }

    Ok(())
}
