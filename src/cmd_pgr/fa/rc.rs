use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};
use std::collections::HashSet;
use std::io::Write;

/// Build the clap subcommand for rc.
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
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input FASTA file to process",
        ))
        .arg(crate::cmd_pgr::args::fa_name_list_arg(false))
        .arg(
            Arg::new("consistent")
                .long("consistent")
                .short('c')
                .action(ArgAction::SetTrue)
                .help("Keep the name consistent (don't prepend 'RC_')"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the rc command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut protected: Vec<&str> = vec![infile.as_str()];
    if let Some(list) = args.get_one::<String>("name_list") {
        protected.push(list.as_str());
    }
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, protected)?;
    let mut reader = pgr::libs::fmt::seq::SeqReader::new(infile)
        .with_context(|| format!("Failed to open reader for {}", infile))?;
    let mut rec = pgr::libs::fmt::seq::SeqRecord::new();

    let is_consistent = args.get_flag("consistent");

    let mut fa_out = pgr::libs::fmt::fa::writer(outfile)
        .with_context(|| format!("Failed to open writer for {}", outfile))?;

    let set_list: HashSet<String> = if args.contains_id("name_list") {
        pgr::libs::io::read_names::<std::collections::HashSet<String>>(
            args.get_one::<String>("name_list").unwrap(),
        )?
    } else {
        HashSet::new()
    };

    while reader.read_record(&mut rec)? {
        let name = String::from_utf8(rec.name().to_vec())?;

        if args.contains_id("name_list") && !set_list.contains(&name) {
            let out =
                pgr::libs::fmt::fa::new_record_with_desc(&name, rec.description(), rec.sequence());
            fa_out.write_record(&out)?;
            continue;
        }

        let new_name = if is_consistent {
            name
        } else {
            format!("RC_{}", name)
        };

        // Reverse complement using the `NT_COMP` lookup table. Standard and
        // IUPAC bases are complemented (case preserved); any other byte
        // (e.g. `-`, `*`) is kept as-is, matching the documented behavior
        // ("Non-IUPAC characters are preserved as-is"). The previous
        // `noodles` `Sequence::complement()` errored on such characters.
        let seq_rc: Vec<u8> = rec
            .sequence()
            .iter()
            .rev()
            .map(|&b| {
                let c = pgr::libs::nt::NT_COMP[b as usize];
                if c == 255 {
                    b
                } else {
                    c
                }
            })
            .collect();
        let record_rc =
            pgr::libs::fmt::fa::new_record_with_desc(&new_name, rec.description(), &seq_rc);
        fa_out.write_record(&record_rc)?;
    }

    fa_out.get_mut().flush()?;

    Ok(())
}
