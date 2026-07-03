use clap::*;
use pgr::libs::loc;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("check")
        .about("Checks genome locations in block FA headers")
        .after_help(
            r###"
Checks genome locations in block FA headers against a chrom.sizes file.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

"###,
        )
        .arg(crate::cmd_pgr::args::genome_arg())
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input block FA file(s) to check"),
        )
        .arg(crate::cmd_pgr::args::fas_name_arg(
            "Check sequences for a specific species",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;
    let opt_genome = args.get_one::<String>("genome").unwrap();
    let opt_name = &args
        .get_one::<String>("name")
        .map(|s| s.as_str())
        .unwrap_or("")
        .to_string();

    //----------------------------
    // Ops
    //----------------------------
    let (mut genome_reader, loc_of) = loc::open_indexed(opt_genome, false)?;

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::reader(infile)?;

        while let Ok(block) = pgr::libs::fmt::fas::next_fas_block(&mut reader) {
            let block_names = block.names;

            // Check if a specific species is requested
            if !opt_name.is_empty() && block_names.contains(opt_name) {
                for entry in &block.entries {
                    let entry_name = entry.range().name();
                    if entry_name == opt_name {
                        let status = pgr::libs::fmt::fas::check_entry_against_ref(
                            entry,
                            &mut genome_reader,
                            &loc_of,
                        )?;
                        writer.write_all(format!("{}\t{}\n", entry.range(), status).as_ref())?;
                    }
                }
            } else if opt_name.is_empty() {
                // Check all sequences in the block
                for entry in &block.entries {
                    let status = pgr::libs::fmt::fas::check_entry_against_ref(
                        entry,
                        &mut genome_reader,
                        &loc_of,
                    )?;
                    writer.write_all(format!("{}\t{}\n", entry.range(), status).as_ref())?;
                }
            }
        }
    }

    Ok(())
}
