use anyhow::Context;
use clap::{ArgMatches, Command};
use std::io::Write;

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
   pgr fas create tests/fas/I.connect.tsv -g tests/fas/genome.fa

2. Create block FA files for a specific species:
   pgr fas create tests/fas/I.connect.tsv -g tests/fas/genome.fa --name S288c

"###,
        )
        .arg(crate::cmd_pgr::args::genome_arg())
        .arg(crate::cmd_pgr::args::infiles_arg_with_help(
            "Input file(s) containing links of ranges",
        ))
        .arg(crate::cmd_pgr::args::fas_name_arg(
            "Set a species name for ranges (default: inferred from header)",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the create command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let opt_genome = args.get_one::<String>("genome").unwrap();
    let opt_name: &str = args
        .get_one::<String>("name")
        .map(|s| s.as_str())
        .unwrap_or("");

    // The writer is opened (truncating) before the input links and the
    // reference genome are read, so reject an `-o` that matches any of them.
    // The genome's `.loc` sidecar is also protected: truncating it would make
    // `open_indexed` (inside `create_from_links`) treat it as fresh and serve
    // an empty index, silently dropping every link.
    let mut inputs: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .map(|s| s.to_string())
        .collect();
    inputs.push(opt_genome.to_string());
    inputs.push(format!("{}.loc", opt_genome));
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, inputs.iter().map(|s| s.as_str()))?;

    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader =
            pgr::reader(infile).with_context(|| format!("Failed to open reader for {}", infile))?;
        pgr::libs::fmt::fas::create_from_links(reader, &mut writer, opt_genome, opt_name)?;
    }

    writer.flush()?;
    Ok(())
}
