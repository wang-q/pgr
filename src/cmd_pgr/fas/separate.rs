use clap::{Arg, ArgAction, ArgMatches, Command};
use std::collections::BTreeMap;
use std::io::Write;

/// Build the clap subcommand for separate.
pub fn make_subcommand() -> Command {
    Command::new("separate")
        .about("Separates block FA files by species")
        .after_help(
            r###"
Separates block FA files by species, creating individual output files for each species.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* Dashes ('-') in sequences are removed
* If the target file already exists, it will be overwritten
* Optionally, sequences can be reverse-complemented if the chromosome strand is '-'

Examples:
1. Separate block FA files by species:
   pgr fas separate tests/fas/example.fas -o output_dir

2. Separate block FA files and reverse-complement sequences:
   pgr fas separate tests/fas/example.fas -o output_dir --rc

3. Use a custom suffix for output files:
   pgr fas separate tests/fas/example.fas -o output_dir --suffix .fa

4. Output to stdout:
   pgr fas separate tests/fas/example.fas

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(crate::cmd_pgr::args::suffix_arg(".fasta"))
        .arg(
            Arg::new("rc")
                .long("rc")
                .action(ArgAction::SetTrue)
                .help("Reverse-complement sequences when chromosome strand is '-'"),
        )
        .arg(crate::cmd_pgr::args::outdir_arg())
}

/// Execute the separate command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outdir = args.get_one::<String>("outdir").unwrap();
    if outdir != "stdout" {
        std::fs::create_dir_all(outdir)?;
    }

    let opt_suffix = args.get_one::<String>("suffix").unwrap();
    let is_rc = args.get_flag("rc");

    let mut file_of: BTreeMap<String, std::fs::File> = BTreeMap::new();
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::reader(infile)?;

        while let Ok(block) = pgr::libs::fmt::fas::next_fas_block(&mut reader) {
            for entry in &block.entries {
                let entry_name = entry.range().name(); // Don't borrow the following `range`
                let mut range = entry.range().clone();

                // Reverse-complement the sequence if needed
                let seq = if is_rc && range.strand() == "-" {
                    *range.strand_mut() = "+".to_string();
                    pgr::libs::nt::rev_comp(entry.seq()).collect::<Vec<u8>>()
                } else {
                    entry.seq().to_vec()
                };

                // Remove dashes from the sequence
                let seq = std::str::from_utf8(&seq)?.to_string().replace('-', "");

                if outdir == "stdout" {
                    print!(">{}\n{}\n", range, seq);
                } else {
                    if !file_of.contains_key(entry_name) {
                        let path =
                            std::path::Path::new(outdir).join(range.name().to_owned() + opt_suffix);
                        let file = std::fs::OpenOptions::new()
                            .create(true)
                            .write(true)
                            .truncate(true)
                            .open(path)?;
                        file_of.insert(entry_name.to_string(), file);
                    }
                    write!(
                        file_of.get(entry_name).ok_or_else(|| anyhow::anyhow!(
                            "file not found for entry: {}",
                            entry_name
                        ))?,
                        ">{}\n{}\n",
                        range,
                        seq
                    )?;
                }
            }
        }
    }

    Ok(())
}
