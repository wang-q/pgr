use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};
use std::collections::BTreeMap;
use std::io::Write;

/// Build the clap subcommand for concat.
pub fn make_subcommand() -> Command {
    Command::new("concat")
        .about("Concatenates sequence pieces of the same species")
        .after_help(
            r###"
Concatenates sequence pieces of the same species from block FA files.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* The --required file lists species names to keep, one per line
* The order of species in the output follows the order in the <name.lst> file
* Missing sequences are filled with gaps (`-`)

Examples:
1. Concatenate sequences and output in FASTA format:
   pgr fas concat tests/fas/example.fas -R tests/fas/name.lst

2. Concatenate sequences and output in relaxed PHYLIP format:
   pgr fas concat tests/fas/example.fas -R tests/fas/name.lst --phylip

3. Output results to a file:
   pgr fas concat tests/fas/example.fas -R tests/fas/name.lst -o output.fas

"###,
        )
        .arg(crate::cmd_pgr::args::required_species_list_arg())
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(
            Arg::new("phylip")
                .long("phylip")
                .action(ArgAction::SetTrue)
                .help("Output in relaxed PHYLIP format instead of FA"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the concat command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
    let is_phylip = args.get_flag("phylip");

    let needed =
        pgr::libs::io::read_names::<Vec<String>>(args.get_one::<String>("required").unwrap())?;
    anyhow::ensure!(!needed.is_empty(), "--required file is empty");

    let mut seq_of: BTreeMap<String, String> = BTreeMap::new();
    for name in &needed {
        seq_of.insert(name.to_string(), String::new());
    }

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader =
            pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;
        pgr::libs::fmt::fas::concat_blocks_into(&mut reader, &needed, &mut seq_of)?;
    }

    pgr::libs::fmt::fas::write_concat_output(&mut writer, &needed, &seq_of, is_phylip)?;
    writer.flush()?;
    Ok(())
}
