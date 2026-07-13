use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};
use std::collections::BTreeMap;
use std::io::{BufWriter, Write};

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

    let mut file_of: BTreeMap<String, BufWriter<std::fs::File>> = BTreeMap::new();
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader =
            pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;

        for block_result in pgr::libs::fmt::fas::iter_fas_blocks(&mut reader) {
            let block = block_result?;
            for (idx, entry) in block.entries.iter().enumerate() {
                let entry_name = &block.names[idx];

                // Reverse-complement the sequence if needed
                let (range_str, seq) = if is_rc && entry.range().strand() == "-" {
                    let mut range = entry.range().clone();
                    *range.strand_mut() = "+".to_string();
                    (
                        range.to_string(),
                        pgr::libs::nt::rev_comp(entry.seq()).collect::<Vec<u8>>(),
                    )
                } else {
                    (entry.range().to_string(), entry.seq().to_vec())
                };

                // Remove dashes from the sequence
                let seq = std::str::from_utf8(&seq)?.replace('-', "");

                if outdir == "stdout" {
                    write!(out, ">{}\n{}\n", range_str, seq)?;
                } else {
                    let file_key = pgr::libs::io::sanitize_filename(entry_name);
                    if !file_of.contains_key(&file_key) {
                        let path = std::path::Path::new(outdir)
                            .join(format!("{}{}", file_key, opt_suffix));
                        let file = std::fs::OpenOptions::new()
                            .create(true)
                            .write(true)
                            .truncate(true)
                            .open(path)?;
                        file_of.insert(file_key.clone(), BufWriter::new(file));
                    }
                    write!(
                        file_of.get_mut(&file_key).unwrap(),
                        ">{}\n{}\n",
                        range_str,
                        seq
                    )?;
                }
            }
        }
    }

    // Explicitly flush all file handles to catch errors on close (e.g. disk full)
    for fh in file_of.values_mut() {
        fh.flush()?;
    }
    out.flush()?;

    Ok(())
}
