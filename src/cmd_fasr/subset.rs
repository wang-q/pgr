use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("subset")
        .about("Extract a subset of species from block FA files")
        .after_help(
            r###"
Extract a subset of species from block FA files based on a list of names.

* <name.lst>: A file containing a list of species names to keep, one per line.
    - The order of species in the output will follow the order in <name.lst>.

* <name.lst> is a file with a list of names to keep, one per line
    * Orders in the output file will following the ones in <name.lst>

* <infiles> are paths to block fasta files, .fas.gz is supported
    * infile == stdin means reading from STDIN

"###,
        )
        .arg(
            Arg::new("name.lst")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Path to name.lst"),
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(2)
                .help("Set the input files to use"),
        )
        .arg(
            Arg::new("required")
                .long("required")
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
    let is_required = args.get_flag("required");

    let needed = intspan::read_first_column(args.get_one::<String>("name.lst").unwrap());

    //----------------------------
    // Operating
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = intspan::reader(infile);

        'BLOCK: while let Ok(block) = pgr::next_fas_block(&mut reader) {
            let block_names = block.names;

            if is_required {
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
