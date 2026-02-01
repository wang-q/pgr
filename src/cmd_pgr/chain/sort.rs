use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::chain::read_chains;
use std::cmp::Ordering;
use std::fs::File;
use std::io::BufWriter;

pub fn make_subcommand() -> Command {
    Command::new("sort")
        .about("Sort chains by score")
        .arg(
            Arg::new("files")
                .required(true)
                .num_args(1..)
                .action(ArgAction::Append)
                .help("Input chain file(s)"),
        )
        .arg(
            Arg::new("output")
                .long("output")
                .short('o')
                .help("Output file (default: stdout)"),
        )
        .arg(
            Arg::new("save_id")
                .long("save-id")
                .action(ArgAction::SetTrue)
                .help("Keep existing chain IDs (default: renumber starting from 1)"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let files: Vec<&String> = args.get_many("files").unwrap().collect();
    let output = args.get_one::<String>("output");
    let save_id = args.get_flag("save_id");

    let mut all_chains = Vec::new();

    // Read all chains
    for file_path in files {
        let file = File::open(file_path)?;
        let chains = read_chains(file)?;
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
    if let Some(out_path) = output {
        let mut writer = BufWriter::new(File::create(out_path)?);
        for chain in all_chains {
            chain.write(&mut writer)?;
        }
    } else {
        let mut writer = BufWriter::new(std::io::stdout());
        for chain in all_chains {
            chain.write(&mut writer)?;
        }
    }

    Ok(())
}
