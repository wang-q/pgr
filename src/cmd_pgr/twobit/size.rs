use anyhow::Context;
use clap::{ArgMatches, Command};
use pgr::libs::fmt::twobit::TwoBitFile;
use std::io::Write;

/// Build the clap subcommand for size.
pub fn make_subcommand() -> Command {
    Command::new("size")
        .about("Counts total bases in 2bit file(s)")
        .after_help(
            r###"
This command retrieves the sequence sizes from a 2bit file.

Notes:
* 2bit files are binary and require random access (seeking)
* Does not support stdin or gzipped inputs

Examples:
1. Get sizes from a 2bit file:
   pgr 2bit size input.2bit

2. Save the output to a file:
   pgr 2bit size input.2bit -o output.tsv

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("2bit"))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(crate::cmd_pgr::args::no_ns_arg())
}

/// Execute the size command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    let no_ns = args.get_flag("no_ns");

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut tb = TwoBitFile::open(infile)
            .with_context(|| format!("Failed to open 2bit file {}", infile))?;

        // Get all sequence names
        let mut names = tb.get_sequence_names();
        // Sort names to be deterministic (optional but good)
        names.sort();

        for name in names {
            if no_ns {
                let len = tb.get_sequence_len_no_ns(&name)?;
                writer.write_fmt(format_args!("{}\t{}\n", name, len))?;
            } else {
                let len = tb.get_sequence_len(&name)?;
                writer.write_fmt(format_args!("{}\t{}\n", name, len))?;
            }
        }
    }

    Ok(())
}
