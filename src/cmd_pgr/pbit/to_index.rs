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
`pgr dist pgi`. With multiple references, `--ref` selects one by name or
index (default 0).

Notes:
* Archives created without --index contain no index segment and cause an error.

Examples:
1. Extract the embedded reference index:
   pgr pbit to-index cohort.pbit -o ref.pgi

2. Extract a specific reference from a multi-reference archive:
   pgr pbit to-index cohort.pbit --ref 1 -o ref2.pgi
"###,
        )
        .arg(
            Arg::new("infile")
                .index(1)
                .required(true)
                .help("Input pbit archive"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg_required())
        .arg(
            Arg::new("ref")
                .long("ref")
                .help("Reference name or index to extract (default 0)"),
        )
}

/// Execute the to-index command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = args.get_one::<String>("outfile").unwrap();
    let ref_spec = args.get_one::<String>("ref").map(|s| s.as_str());
    let mut dec = pgr::libs::pbit::decompressor::Decompressor::open(infile)?;
    let ref_names: Vec<String> = dec.ref_table().iter().map(|r| r.ref_name.clone()).collect();
    let ref_id = match ref_spec {
        Some(s) => {
            if let Ok(n) = s.parse::<usize>() {
                anyhow::ensure!(
                    n < ref_names.len(),
                    "reference index {} out of range ({} references)",
                    n,
                    ref_names.len()
                );
                n
            } else {
                ref_names.iter().position(|n| n == s).ok_or_else(|| {
                    anyhow::anyhow!("reference '{}' not found (available: {ref_names:?})", s)
                })?
            }
        }
        None => 0,
    };
    let idx = dec
        .read_reference_index(ref_id)?
        .ok_or_else(|| anyhow::anyhow!("no embedded reference index in {}", infile))?;
    let mut writer = pgr::writer(outfile)?;
    idx.write(&mut writer)?;
    log::info!("wrote embedded reference index to {}", outfile);
    Ok(())
}
