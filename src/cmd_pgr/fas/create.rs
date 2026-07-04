use clap::{ArgMatches, Command};

/// Build the clap subcommand for create.
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

/// Execute the create command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;
    let opt_genome = args.get_one::<String>("genome").unwrap();
    let opt_name = &args
        .get_one::<String>("name")
        .map(|s| s.as_str())
        .unwrap_or("")
        .to_string();

    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = pgr::reader(infile)?;
        pgr::libs::fmt::fas::create_from_links(reader, &mut writer, opt_genome, opt_name)?;
    }

    Ok(())
}
