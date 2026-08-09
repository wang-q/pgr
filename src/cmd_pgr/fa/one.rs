use anyhow::Context;
use clap::{Arg, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for one.
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
            Arg::new("seq_name")
                .required(true)
                .index(2)
                .help("Name of the sequence to extract"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the one command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, [infile.as_str()])?;
    let mut reader = pgr::libs::fmt::seq::SeqReader::new(infile)
        .with_context(|| format!("Failed to open reader for {}", infile))?;
    let mut rec = pgr::libs::fmt::seq::SeqRecord::new();

    let mut fa_out = pgr::libs::fmt::fa::writer(outfile)
        .with_context(|| format!("Failed to open writer for {}", outfile))?;

    let name = args.get_one::<String>("seq_name").unwrap();

    let mut found = false;
    while reader.read_record(&mut rec)? {
        if rec.name() == name.as_bytes() {
            let out = pgr::libs::fmt::fa::new_record_with_desc(
                &String::from_utf8_lossy(rec.name()),
                rec.description(),
                rec.sequence(),
            );
            fa_out.write_record(&out)?;
            found = true;
            break;
        }
    }

    if !found {
        anyhow::bail!("sequence {} not found in {}", name, infile);
    }

    fa_out.get_mut().flush()?;

    Ok(())
}
