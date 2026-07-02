use clap::*;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("subset")
        .about("Extracts a subset of species from block FA files")
        .after_help(
            r###"
Extracts a subset of species from block FA files.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* The --required file lists species names to keep, one per line
* The order of species in the output follows the order in the <name.lst> file

Examples:
1. Extract a subset of species:
   pgr fas subset tests/fas/example.fas -R tests/fas/name.lst

2. Extract a subset and skip blocks missing any required species:
   pgr fas subset tests/fas/example.fas -R tests/fas/name.lst --strict

3. Output results to a file:
   pgr fas subset tests/fas/example.fas -R tests/fas/name.lst -o output.fas

"###,
        )
        .arg(
            Arg::new("required")
                .long("required")
                .short('R')
                .required(true)
                .num_args(1)
                .help("File with a list of species names to keep, one per line"),
        )
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(
            Arg::new("strict")
                .long("strict")
                .action(ArgAction::SetTrue)
                .help("Skip blocks not containing all the names"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;
    let is_strict = args.get_flag("strict");

    let needed =
        pgr::libs::io::read_names::<Vec<String>>(args.get_one::<String>("required").unwrap())?;

    //----------------------------
    // Operating
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::reader(infile)?;

        'BLOCK: while let Ok(block) = pgr::libs::fmt::fas::next_fas_block(&mut reader) {
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
