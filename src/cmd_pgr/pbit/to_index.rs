//! Extract the embedded reference index (.pgi) from a pbit archive.

use clap::{Arg, ArgMatches, Command};

/// Build the clap subcommand for to-index.
pub fn make_subcommand() -> Command {
    Command::new("to-index")
        .about("Extracts the embedded reference .pgi index from a pbit archive")
        .after_help(
            r###"
Reads the reference index segment embedded by `pgr pbit create --index` and
writes it as a standalone `.pgi` file, usable by `pgr pgi align` /
`pgr dist pgi`.

Notes:
* Archives created without --index contain no index segment and cause an error.

Examples:
1. Extract the embedded reference index:
   pgr pbit to-index cohort.pbit -o ref.pgi
"###,
        )
        .arg(
            Arg::new("infile")
                .index(1)
                .required(true)
                .help("Input pbit archive"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg_required())
}

/// Execute the to-index command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = args.get_one::<String>("outfile").unwrap();
    let mut dec = pgr::libs::pbit::decompressor::Decompressor::open(infile)?;
    let idx = dec
        .read_reference_index()?
        .ok_or_else(|| anyhow::anyhow!("no embedded reference index in {}", infile))?;
    let mut writer = pgr::writer(outfile)?;
    idx.write(&mut writer)?;
    log::info!("wrote embedded reference index to {}", outfile);
    Ok(())
}
