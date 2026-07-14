use anyhow::Context;
use clap::{Arg, ArgMatches, Command};
use pgr::libs::chain::net::{write_net_file, ChainNet};
use pgr::libs::chain::ChainReader;

/// Build the clap subcommand for net.
pub fn make_subcommand() -> Command {
    Command::new("net")
        .about("Makes alignment nets out of chains")
        .after_help(
            r###"
Builds alignment nets from chains. Nets hierarchically organize chains into
filled alignments and gaps, providing a layered view of synteny between two
genomes.

Notes:
* Input chain file must already be sorted by score descending (use `pgr chain sort`); otherwise the command returns an error
* Outputs two net files: one in target orientation, one in query orientation
* Use `--min-space` to control the minimum gap size to fill (default: 25)
* Use `--min-fill` to control the minimum fill to record (default: min-space / 2)
* Use `--min-score` to filter low-scoring chains (default: 2000)
* Use `--incl-hap` to include haplotype chains (names containing `_hap` or `_alt`)

Examples:
1. Build nets from sorted chains:
   pgr chain net in.chain t.sizes q.sizes t.net q.net

2. Adjust fill parameters:
   pgr chain net in.chain t.sizes q.sizes t.net q.net --min-space 50 --min-fill 20

3. Include haplotype chains:
   pgr chain net in.chain t.sizes q.sizes t.net q.net --incl-hap

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input chain file",
        ))
        .arg(crate::cmd_pgr::args::chain_t_sizes_arg())
        .arg(crate::cmd_pgr::args::chain_q_sizes_arg())
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
        .arg(crate::cmd_pgr::args::incl_hap_arg())
}
/// Execute the net command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
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

    let mut reader = ChainReader::new(
        pgr::reader(input_path)
            .with_context(|| format!("Failed to open reader for {}", input_path))?,
    );

    let mut last_score = f64::MAX;

    for res in reader.by_ref() {
        let chain = res?;

        // Input must be sorted by score descending.
        if chain.header.score > last_score {
            anyhow::bail!(
                "Input not sorted by score: {} > {}",
                chain.header.score,
                last_score
            );
        }
        last_score = chain.header.score;

        if chain.header.score < min_score {
            continue;
        }

        if !incl_hap && pgr::libs::chain::pre_net::is_haplotype(&chain.header.q_name) {
            continue;
        }

        // Add to T net
        t_net.add_chain(chain.clone(), min_space, min_fill, min_score);

        // Add to Q net
        q_net.add_chain_as_q(chain, min_space, min_fill, min_score);
    }

    // Finish and write T net
    write_net_file(
        target_net_path,
        &t_net,
        false,
        &reader.header_comments,
        min_score,
        min_fill,
    )?;

    // Finish and write Q net
    write_net_file(
        query_net_path,
        &q_net,
        true,
        &reader.header_comments,
        min_score,
        min_fill,
    )?;

    Ok(())
}
