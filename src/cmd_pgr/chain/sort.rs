use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::chain::read_chains;
use std::cmp::Ordering;
use std::io::BufRead;
/// Build the clap subcommand for sort.
pub fn make_subcommand() -> Command {
    Command::new("sort")
        .about("Sorts chains by score")
        .arg(
            Arg::new("infiles")
                .required_unless_present("input_list")
                .num_args(1..)
                .action(ArgAction::Append)
                .help("Input chain file(s)"),
        )
        .arg(
            Arg::new("input_list")
                .long("input-list")
                .help("File containing a list of input chain files (one per line)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("save_id")
                .long("save-id")
                .action(ArgAction::SetTrue)
                .help("Keep existing chain IDs (default: renumber starting from 1)"),
        )
}
/// Execute the sort command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut files: Vec<String> = args
        .get_many::<String>("infiles")
        .map(|v| v.cloned().collect())
        .unwrap_or_default();

    if let Some(list_path) = args.get_one::<String>("input_list") {
        let reader = pgr::reader(list_path)?;
        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if !trimmed.is_empty() {
                files.push(trimmed.to_string());
            }
        }
    }

    let save_id = args.get_flag("save_id");

    let mut all_chains = Vec::new();

    // Read all chains
    for file_path in &files {
        let chains = read_chains(pgr::reader(file_path)?)?;
        all_chains.extend(chains);
    }

    // Sort by score descending
    all_chains.sort_by(|a, b| {
        b.header
            .score
            .partial_cmp(&a.header.score)
            .unwrap_or(Ordering::Equal)
    });

    // Renumber if needed
    if !save_id {
        for (i, chain) in all_chains.iter_mut().enumerate() {
            chain.header.id = (i + 1) as u64;
        }
    }

    // Write output
    let out_path = crate::cmd_pgr::args::get_outfile(args);
    let mut writer = pgr::writer(out_path)?;
    for chain in all_chains {
        chain.write(&mut writer)?;
    }

    Ok(())
}
