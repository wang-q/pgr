use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for replace.
pub fn make_subcommand() -> Command {
    Command::new("replace")
        .about("Replaces headers of a FASTA file based on a TSV mapping")
        .after_help(
            r###"
This command replaces sequence headers in a FASTA file based on a TSV mapping file.
The TSV file format:
    seq1    replace_name    more_replace_name
    seq2    replace_name
    seq2    another_replace_name

Notes:
* The TSV file should contain two or more columns: the original name and the replacement name
* If more than two columns are provided, the sequence will be duplicated for each replacement name
* Multiple lines of the same original_name will also duplicate the record
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

Examples:
1. Replace headers using a TSV file:
   pgr fa replace input.fa --replace-tsv replace.tsv -o output.fa

2. Only output sequences listed in the TSV file (like `pgr fa some`):
   pgr fa replace input.fa --replace-tsv replace.tsv --some -o output.fa

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input FASTA file to process",
        ))
        .arg(crate::cmd_pgr::args::replace_tsv_arg())
        .arg(
            Arg::new("some")
                .long("some")
                .action(ArgAction::SetTrue)
                .help("Only output sequences listed in the TSV file, like `pgr fa some`"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the replace command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let mut fa_in = pgr::libs::fmt::fa::reader(infile)
        .with_context(|| format!("Failed to open reader for {}", infile))?;

    let replace_of =
        pgr::libs::io::read_replace_tsv(args.get_one::<String>("replace_tsv").unwrap())?;
    let is_some = args.get_flag("some");

    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut fa_out = pgr::libs::fmt::fa::writer(outfile)
        .with_context(|| format!("Failed to open writer for {}", outfile))?;

    for result in fa_in.records() {
        let record = result?;
        let name = String::from_utf8(record.name().into())?;

        if let Some(new_names) = replace_of.get(&name) {
            for el in new_names {
                let record_replace = pgr::libs::fmt::fa::new_record_preserving_desc(
                    el,
                    &record,
                    &record.sequence()[..],
                );
                fa_out.write_record(&record_replace)?;
            }
        } else if !is_some {
            fa_out.write_record(&record)?;
        }
    }

    fa_out.get_mut().flush()?;

    Ok(())
}
