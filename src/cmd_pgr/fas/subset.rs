use clap::*;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("subset")
        .about("Extract a subset of species from block FA files")
        .after_help(
            r###"
Extracts a subset of species from block FA files.

Input files can be gzipped. If the input file is 'stdin', data is read from standard input.

Note:
- The --required file lists species names to keep, one per line.
- The order of species in the output follows the order in the <name.lst> file.

Examples:
1. Extract a subset of species:
   pgr fas subset tests/fas/example.fas -r tests/fas/name.lst

2. Extract a subset and skip blocks missing any required species:
   pgr fas subset tests/fas/example.fas -r tests/fas/name.lst --strict

3. Output results to a file:
   pgr fas subset tests/fas/example.fas -r tests/fas/name.lst -o output.fas

"###,
        )
        .arg(
            Arg::new("name.lst")
                .short('r')
                .long("required")
                .required(true)
                .num_args(1)
                .help("Required: File with a list of species names to keep"),
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input block FA file(s) to process"),
        )
        .arg(
            Arg::new("strict")
                .long("strict")
                .action(ArgAction::SetTrue)
                .help("Skip blocks not containing all the names"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let is_strict = args.get_flag("strict");

    let needed = intspan::read_first_column(args.get_one::<String>("name.lst").unwrap());

    //----------------------------
    // Operating
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = intspan::reader(infile);

        'BLOCK: while let Ok(block) = pgr::libs::fas::next_fas_block(&mut reader) {
            let block_names = block.names;

            if is_strict {
                for name in &needed {
                    if !block_names.contains(name) {
                        continue 'BLOCK;
                    }
                }
            }

            for name in &needed {
                if block_names.contains(name) {
                    for entry in &block.entries {
                        let entry_name = entry.range().name();
                        //----------------------------
                        // Output
                        //----------------------------
                        if entry_name == name {
                            writer.write_all(entry.to_string().as_ref())?;
                        }
                    }
                }
            }

            // end of a block
            writer.write_all("\n".as_ref())?;
        }
    }

    Ok(())
}
