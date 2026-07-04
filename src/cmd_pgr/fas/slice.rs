use clap::{ArgMatches, Command};

/// Build the clap subcommand for slice.
pub fn make_subcommand() -> Command {
    Command::new("slice")
        .about("Extracts alignment slices")
        .after_help(
            r###"
Extracts alignment slices from block FA files using a runlist JSON.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* The JSON file (--runlist) keys are chromosome/sequence names, and values are runlists (e.g., "1-100,200-300")

Examples:
1. Extract slices defined in a JSON file:
   pgr fas slice tests/fas/slice.fas --runlist tests/fas/slice.json

2. Extract slices and name the output based on a specific species:
   pgr fas slice tests/fas/slice.fas --runlist tests/fas/slice.json --name S288c

3. Output results to a file:
   pgr fas slice tests/fas/slice.fas --runlist tests/fas/slice.json -o output.fas

"###,
        )
        .arg(crate::cmd_pgr::args::runlist_arg())
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(
            crate::cmd_pgr::args::fas_name_arg("Reference species name. Default is the first species"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the slice command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    let json = intspan::read_json(args.get_one::<String>("runlist").unwrap());
    let set = intspan::json2set(&json);

    let mut name = args
        .get_one::<String>("name")
        .map(|s| s.to_string())
        .unwrap_or_default();

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::reader(infile)?;

        loop {
            let block = match pgr::libs::fmt::fas::next_fas_block(&mut reader) {
                Ok(b) => b,
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(e.into()),
            };
            // the first name of the first block becomes the default reference
            if name.is_empty() {
                name = block.names.first().cloned().unwrap_or_default();
            }

            pgr::libs::alignment::slice_block(&block, &name, &set, &mut writer)?;
        }
    }

    Ok(())
}
