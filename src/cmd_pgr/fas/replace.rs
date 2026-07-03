use clap::{ArgMatches, Command};
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
* The replacement file (--replace-tsv) should contain one or more fields:
  * `original_name  replace_name   more_replace_name`
* One field: Deletes the entire alignment block for the specified species
* Three or more fields: Duplicates the entire alignment block for each replacement name

Examples:
1. Replace species names in a block FA file:
   pgr fas replace tests/fas/example.fas --replace-tsv tests/fas/replace.tsv

2. Output results to a file:
   pgr fas replace tests/fas/example.fas --replace-tsv tests/fas/replace.tsv -o output.fas

"###,
        )
        .arg(crate::cmd_pgr::args::replace_tsv_arg())
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;

    let replace_of =
        pgr::libs::io::read_replace_tsv(args.get_one::<String>("replace_tsv").unwrap())?;

    //----------------------------
    // Operating
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::reader(infile)?;

        while let Ok(block) = pgr::libs::fmt::fas::next_fas_block(&mut reader) {
            let blocks = pgr::libs::fmt::fas::replace_block_lines(&block, &replace_of)?;
            for b in &blocks {
                writer.write_all(b.as_ref())?;
                writer.write_all("\n".as_ref())?;
            }
        }
    }

    Ok(())
}
