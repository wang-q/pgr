use clap::*;
use std::collections::BTreeMap;
use std::io::Write;
// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("concat")
        .about("Concatenates sequence pieces of the same species")
        .after_help(
            r###"
Concatenates sequence pieces of the same species from block FA files.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* The --required file lists species names to keep, one per line
* The order of species in the output follows the order in the <name.lst> file
* Missing sequences are filled with gaps (`-`)

Examples:
1. Concatenate sequences and output in FASTA format:
   pgr fas concat tests/fasr/example.fas -r tests/fasr/name.lst

2. Concatenate sequences and output in relaxed PHYLIP format:
   pgr fas concat tests/fasr/example.fas -r tests/fasr/name.lst --phylip

3. Output results to a file:
   pgr fas concat tests/fasr/example.fas -r tests/fasr/name.lst -o output.fas

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
        .arg(
            Arg::new("phylip")
                .long("phylip")
                .action(ArgAction::SetTrue)
                .help("Output in relaxed PHYLIP format instead of FA"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap())?;
    let is_phylip = args.get_flag("phylip");

    let needed = pgr::libs::io::read_names_as_vec(args.get_one::<String>("required").unwrap())?;

    let mut seq_of: BTreeMap<String, String> = BTreeMap::new();
    for name in &needed {
        // default value
        seq_of.insert(name.to_string(), String::new());
    }

    //----------------------------
    // Ops
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::reader(infile)?;

        while let Ok(block) = pgr::libs::fmt::fas::next_fas_block(&mut reader) {
            let block_names = block.names;
            let length = block.entries.first().unwrap().seq().len();

            for name in &needed {
                if block_names.contains(name) {
                    for entry in &block.entries {
                        let entry_name = entry.range().name();
                        if entry_name == name {
                            let seq = std::str::from_utf8(entry.seq()).unwrap();
                            seq_of.entry(name.to_string()).and_modify(|e| *e += seq);
                        }
                    }
                } else {
                    // fill absent names with ------
                    seq_of
                        .entry(name.to_string())
                        .and_modify(|e| *e += "-".repeat(length).as_str());
                }
            }
        }
    }

    //----------------------------
    // Output
    //----------------------------
    if is_phylip {
        let count = needed.len();
        let length = seq_of.first_key_value().unwrap().1.len();
        writer.write_all(format!("{} {}\n", count, length).as_ref())?;
        for (k, v) in &seq_of {
            writer.write_all(format!("{} {}\n", k, v).as_ref())?;
        }
    } else {
        for (k, v) in &seq_of {
            writer.write_all(format!(">{}\n{}\n", k, v).as_ref())?;
        }
    }

    Ok(())
}
