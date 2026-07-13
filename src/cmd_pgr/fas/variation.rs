use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for variation.
pub fn make_subcommand() -> Command {
    Command::new("variation")
        .about("Lists variations (substitutions)")
        .after_help(
            r###"
Lists variations (substitutions) from block FA files in TSV format.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* `--outgroup` requires at least 2 sequences per block and polarizes substitutions against the last sequence
* Filter out complex variations: `tsv-filter -H --ne freq:-1`
* Filter out singletons: `tsv-filter -H --ne freq:1`

Examples:
1. List substitutions from block FA files:
   pgr fas variation tests/fas/example.fas

2. Handle outgroup (last sequence) for polarization:
   pgr fas variation tests/fas/example.fas --outgroup

3. Output results to a file:
   pgr fas variation tests/fas/example.fas -o output.tsv

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(crate::cmd_pgr::args::outgroup_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the variation command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    let has_outgroup = args.get_flag("outgroup");

    let field_names = [
        "#target",
        "chr",
        "chr_pos",
        "range",
        "pos",
        "tbase",
        "qbase",
        "bases",
        "mutant_to",
        "freq",
        "pattern",
        "obase",
    ];

    writeln!(writer, "{}", field_names.join("\t"))?;

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader =
            pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;

        for block_result in pgr::libs::fmt::fas::iter_fas_blocks(&mut reader) {
            let block = block_result?;
            pgr::libs::fmt::fas::write_variations(&block, has_outgroup, &mut writer)?;
        }
    }

    writer.flush()?;
    Ok(())
}
