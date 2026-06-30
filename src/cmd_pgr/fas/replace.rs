use clap::*;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("replace")
        .about("Replaces headers in block FA files")
        .after_help(
            r###"
Replaces headers in block FA files using a TSV file.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* The replacement file (--required) should contain one or more fields:
  * `original_name  replace_name   more_replace_name`
* One field: Deletes the entire alignment block for the specified species
* Three or more fields: Duplicates the entire alignment block for each replacement name

Examples:
1. Replace species names in a block FA file:
   pgr fas replace tests/fas/example.fas -r tests/fas/replace.tsv

2. Output results to a file:
   pgr fas replace tests/fas/example.fas -r tests/fas/replace.tsv -o output.fas

"###,
        )
        .arg(crate::cmd_pgr::fas::common::required_arg())
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input block FA file(s) to process"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap())?;

    let replace_of = pgr::libs::io::read_replace_tsv(args.get_one::<String>("required").unwrap())?;

    //----------------------------
    // Operating
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::reader(infile)?;

        while let Ok(block) = pgr::libs::fmt::fas::next_fas_block(&mut reader) {
            let originals = block.headers.clone();

            let matched: Vec<String> = replace_of
                .keys()
                .filter(|e| originals.contains(*e))
                .map(|e| e.to_string())
                .collect();

            if matched.is_empty() || matched.len() > 1 {
                // block untouched

                if matched.len() > 1 {
                    log::warn!("Doesn't support replacing multiple records in one block");
                }

                //----------------------------
                // Output
                //----------------------------
                for entry in &block.entries {
                    writer.write_all(entry.to_string().as_ref())?;
                }
                writer.write_all("\n".as_ref())?;
            } else {
                let original = matched.first().unwrap();
                let idx = block.headers.iter().position(|e| e == original).unwrap();
                for new in replace_of.get(original).unwrap() {
                    for (i, entry) in block.entries.iter().enumerate() {
                        //----------------------------
                        // Output
                        //----------------------------
                        if i == idx {
                            writer.write_all(
                                format!(
                                    ">{}\n{}\n",
                                    new,
                                    String::from_utf8(entry.seq().to_vec()).unwrap()
                                )
                                .as_ref(),
                            )?;
                        } else {
                            writer.write_all(entry.to_string().as_ref())?;
                        }
                    }

                    writer.write_all("\n".as_ref())?;
                }
            }
        }
    }

    Ok(())
}
