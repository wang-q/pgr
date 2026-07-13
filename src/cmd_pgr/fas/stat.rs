use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for stat.
pub fn make_subcommand() -> Command {
    Command::new("stat")
        .about("Calculates basic statistics of block FA files")
        .after_help(
            r###"
Calculates basic statistics of block FA files (length, comparable, difference, etc.).

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* Output columns: target length comparable difference gap ambiguous D indel
* `--outgroup` excludes the last sequence from all calculations except length

Examples:
1. Get statistics for block FA files:
   pgr fas stat tests/fas/example.fas

2. Statistics treating the last sequence as an outgroup:
   pgr fas stat tests/fas/example.fas --outgroup

3. Output results to a file:
   pgr fas stat tests/fas/example.fas -o output.tsv

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(crate::cmd_pgr::args::outgroup_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the stat command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    let has_outgroup = args.get_flag("outgroup");

    let field_names = [
        "target",
        "length",
        "comparable",
        "difference",
        "gap",
        "ambiguous",
        "D",
        "indel",
    ];

    writeln!(writer, "{}", field_names.join("\t"))?;

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader =
            pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;

        for block_result in pgr::libs::fmt::fas::iter_fas_blocks(&mut reader) {
            let block = block_result?;
            let stat = pgr::libs::fmt::fas::compute_block_stat(&block, has_outgroup)?;
            writeln!(writer, "{}", stat.to_tsv())?;
        }
    }

    writer.flush()?;
    Ok(())
}
