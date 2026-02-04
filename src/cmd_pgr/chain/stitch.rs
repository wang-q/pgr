use anyhow::Result;
use clap::{Arg, Command};
use pgr::libs::chain::{Chain, ChainReader};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};

pub fn make_subcommand() -> Command {
    Command::new("stitch")
        .about("Join chain fragments with the same chain ID into a single chain per ID")
        .arg(Arg::new("input").required(true).help("Input chain file"))
        .arg(Arg::new("output").required(true).help("Output chain file"))
}

pub fn execute(args: &clap::ArgMatches) -> Result<()> {
    let input_path = args.get_one::<String>("input").unwrap();
    let output_path = args.get_one::<String>("output").unwrap();

    let f = File::open(input_path)?;
    let mut reader = ChainReader::new(BufReader::new(f));

    // Store chains by ID
    let mut chains: HashMap<u64, Chain> = HashMap::new();

    while let Some(res) = reader.next() {
        let chain = res?;

        chains
            .entry(chain.header.id)
            .and_modify(|existing| {
                // Merge logic: tackOnFrag
                // Check consistency
                if existing.header.t_name != chain.header.t_name
                    || existing.header.q_name != chain.header.q_name
                    || existing.header.q_strand != chain.header.q_strand
                {
                    eprintln!(
                        "Warning: Inconsistent chain info for ID {}: skipping fragment",
                        chain.header.id
                    );
                    return;
                }

                // Convert both to blocks
                let mut blocks = existing.to_blocks();
                let new_blocks = chain.to_blocks();

                // Append new blocks
                blocks.extend(new_blocks);

                // Sort blocks by t_start (and q_start for stability if needed)
                blocks.sort_by(|a, b| a.t_start.cmp(&b.t_start));

                // Reconstruct data from blocks
                // This updates header ranges automatically
                existing.data = Chain::from_blocks(&mut existing.header, &blocks);

                // Sum score
                existing.header.score += chain.header.score;
            })
            .or_insert(chain);
    }

    // Collect and sort by score (descending)
    let mut chain_list: Vec<Chain> = chains.into_values().collect();
    chain_list.sort_by(|a, b| {
        // Sort by score descending
        b.header.score.partial_cmp(&a.header.score).unwrap()
    });

    // Write output
    let out_file = File::create(output_path)?;
    let mut writer = BufWriter::new(out_file);

    for chain in chain_list {
        chain.write(&mut writer)?;
    }

    Ok(())
}
