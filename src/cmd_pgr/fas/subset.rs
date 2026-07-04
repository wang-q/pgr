use clap::{Arg, ArgAction, ArgMatches, Command};
use std::collections::HashMap;
use std::io::Write;

/// Build the clap subcommand for subset.
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
        .arg(crate::cmd_pgr::args::required_species_list_arg())
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(
            Arg::new("strict")
                .long("strict")
                .action(ArgAction::SetTrue)
                .help("Skip blocks not containing all the names"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the subset command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;
    let is_strict = args.get_flag("strict");

    let needed =
        pgr::libs::io::read_names::<Vec<String>>(args.get_one::<String>("required").unwrap())?;

    // Operating
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::reader(infile)?;

        for block_result in pgr::libs::fmt::fas::iter_fas_blocks(&mut reader) {
            let block = block_result?;

            // Build name -> entry index for O(1) lookup (avoids O(N*M) triple loop)
            let entry_of: HashMap<&str, &pgr::libs::fmt::fas::FasEntry> = block
                .entries
                .iter()
                .map(|e| (e.range().name().as_str(), e))
                .collect();

            if is_strict && !needed.iter().all(|n| entry_of.contains_key(n.as_str())) {
                continue;
            }

            for name in &needed {
                if let Some(entry) = entry_of.get(name.as_str()) {
                    writer.write_all(entry.to_string().as_ref())?;
                }
            }

            // end of a block
            writer.write_all("\n".as_ref())?;
        }
    }

    Ok(())
}
