use anyhow::Result;
use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::chain::net::{finalize_net, write_net, ChainNet};
use pgr::libs::chain::ChainReader;
use std::io::Write;

pub fn make_subcommand() -> Command {
    Command::new("net")
        .about("Make alignment nets out of chains")
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input chain file",
        ))
        .arg(
            Arg::new("t_sizes")
                .required(true)
                .help("Target sequence sizes"),
        )
        .arg(
            Arg::new("q_sizes")
                .required(true)
                .help("Query sequence sizes"),
        )
        .arg(
            Arg::new("out_target_net")
                .required(true)
                .help("Output target net file"),
        )
        .arg(
            Arg::new("out_query_net")
                .required(true)
                .help("Output query net file"),
        )
        .arg(
            Arg::new("min_space")
                .long("min-space")
                .default_value("25")
                .value_parser(clap::value_parser!(u64))
                .help("Minimum gap size to fill"),
        )
        .arg(
            Arg::new("min_fill")
                .long("min-fill")
                .value_parser(clap::value_parser!(u64))
                .help("Minimum fill to record (default: min_space / 2)"),
        )
        .arg(crate::cmd_pgr::args::min_score_arg("2000"))
        .arg(
            Arg::new("incl_hap")
                .long("incl-hap")
                .action(ArgAction::SetTrue)
                .help("Include haplotype pseudochromosome queries"),
        )
}

pub fn execute(args: &ArgMatches) -> Result<()> {
    let input_path = args.get_one::<String>("infile").unwrap();
    let target_sizes_path = args.get_one::<String>("t_sizes").unwrap();
    let query_sizes_path = args.get_one::<String>("q_sizes").unwrap();
    let target_net_path = args.get_one::<String>("out_target_net").unwrap();
    let query_net_path = args.get_one::<String>("out_query_net").unwrap();

    let min_space = *args.get_one::<u64>("min_space").unwrap();
    let min_fill = args
        .get_one::<u64>("min_fill")
        .copied()
        .unwrap_or(min_space / 2);
    let min_score = *args.get_one::<f64>("min_score").unwrap();
    let incl_hap = args.get_flag("incl_hap");

    let t_sizes = pgr::read_sizes::<u64>(target_sizes_path)?;
    let q_sizes = pgr::read_sizes::<u64>(query_sizes_path)?;

    let mut t_net = ChainNet::new(&t_sizes);
    let mut q_net = ChainNet::new(&q_sizes);

    let mut reader = ChainReader::new(pgr::reader(input_path)?);

    let mut last_score = f64::MAX;

    for res in reader.by_ref() {
        let chain = res?;

        // Sort check (optional but good)
        if chain.header.score > last_score {
            // In C code, it doesn't strictly abort, but expects sorted.
            // We can just warn or proceed. The greedy algorithm relies on score sorting.
            // Let's bail if strict, or just log.
            // bail!("Input not sorted by score");
        }
        last_score = chain.header.score;

        if chain.header.score < min_score {
            continue;
        }

        if !incl_hap && is_haplotype(&chain.header.q_name) {
            continue;
        }

        // Add to T net
        t_net.add_chain(chain.clone(), min_space, min_fill, min_score);

        // Add to Q net
        q_net.add_chain_as_q(chain, min_space, min_fill, min_score);
    }

    // Finish and write T net
    {
        let mut writer = pgr::writer(target_net_path)?;

        for comment in &reader.header_comments {
            write!(writer, "{}", comment)?;
        }

        // We need to iterate chroms in order? C code iterates chromList (which is reversed from creation? No, preserved order).
        // Hash map iteration is random.
        // We should sort keys or iterate if we had a list.
        // For deterministic output, let's sort by name.
        let mut t_chrom_names: Vec<_> = t_net.chroms.keys().cloned().collect();
        t_chrom_names.sort();

        for name in t_chrom_names {
            if let Some(chrom_cell) = t_net.chroms.get(&name) {
                let mut chrom = chrom_cell.borrow_mut();
                finalize_net(&mut chrom, false); // is_q = false
                write_net(&chrom, &mut writer, false, min_score, min_fill)?;
            }
        }
    }

    // Finish and write Q net
    {
        let mut writer = pgr::writer(query_net_path)?;

        for comment in &reader.header_comments {
            write!(writer, "{}", comment)?;
        }

        let mut q_chrom_names: Vec<_> = q_net.chroms.keys().cloned().collect();
        q_chrom_names.sort();

        for name in q_chrom_names {
            if let Some(chrom_cell) = q_net.chroms.get(&name) {
                let mut chrom = chrom_cell.borrow_mut();
                finalize_net(&mut chrom, true); // is_q = true
                write_net(&chrom, &mut writer, true, min_score, min_fill)?;
            }
        }
    }

    Ok(())
}

fn is_haplotype(name: &str) -> bool {
    name.contains("_hap") || name.contains("_alt")
}
