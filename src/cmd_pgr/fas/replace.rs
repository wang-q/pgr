use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for replace.
pub fn make_subcommand() -> Command {
    Command::new("replace")
        .about("Replaces headers in block FA files")
        .after_help(
            r###"
Replaces headers in block FA files using a TSV file.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* The replacement file (--replace-tsv) contains tab-separated fields per line:
* One field: if the name uniquely matches one header in a block, the whole block is dropped
* Two fields: `original_name<TAB>new_name` replaces the matching header
* Three or more fields: duplicates the entire alignment block once for every replacement name after the first
* If a block contains multiple matching headers, the block is kept unchanged and a warning is emitted

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

/// Execute the replace command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    let replace_of =
        pgr::libs::io::read_replace_tsv(args.get_one::<String>("replace_tsv").unwrap())?;

    // Operating
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader =
            pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;

        for block_result in pgr::libs::fmt::fas::iter_fas_blocks(&mut reader) {
            let block = block_result?;
            let blocks = pgr::libs::fmt::fas::replace_block_lines(&block, &replace_of)?;
            for b in &blocks {
                writer.write_all(b.as_ref())?;
                writer.write_all("\n".as_ref())?;
            }
        }
    }

    writer.flush()?;
    Ok(())
}
