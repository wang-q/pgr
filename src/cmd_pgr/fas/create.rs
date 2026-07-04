use clap::{ArgMatches, Command};
use std::io::BufRead;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("create")
        .about("Creates block FA files from links of ranges")
        .after_help(
            r###"
Creates block FA files from links of ranges.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* The reference genome(s) must be provided as a multi-sequence FA file, can be bgzipped
* Two styles of FA headers are supported:
  * `>chr` for single-genome self-alignments
  * `>name.chr` for multiple genomes

Examples:
1. Create block FA files for a single genome:
   pgr fas create tests/fasr/I.connect.tsv -g tests/fasr/genome.fa

2. Create block FA files for a specific species:
   pgr fas create tests/fasr/I.connect.tsv -g tests/fasr/genome.fa --name S288c

"###,
        )
        .arg(crate::cmd_pgr::args::genome_arg())
        .arg(crate::cmd_pgr::args::infiles_arg_with_help(
            "Input file(s) containing links of ranges",
        ))
        .arg(crate::cmd_pgr::args::fas_name_arg(
            "Set a species name for ranges. No effects if --multi",
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
    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = pgr::reader(infile)?;
        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split('\t').collect();

            for part in &parts {
                let mut range = intspan::Range::from_str(part);
                if !range.is_valid() {
                    continue;
                }

                // Set the species name if provided
                if !opt_name.is_empty() {
                    *range.name_mut() = opt_name.to_string();
                }

                // Fetch the sequence from the reference genome
                let seq = pgr::libs::loc::get_seq_loc(opt_genome, &range.to_string())?;

                //----------------------------
                // Output
                //----------------------------
                writer.write_all(format!(">{}\n{}\n", range, seq).as_ref())?;
            }

            // End of a block
            writer.write_all("\n".as_ref())?;
        }
    }

    Ok(())
}
