use anyhow::Result;
use clap::{Arg, ArgMatches, Command};
use pgr::libs::chain::BitMap;
use std::collections::HashMap;
/// Build the clap subcommand for pre-net.
pub fn make_subcommand() -> Command {
    Command::new("pre-net")
        .about("Removes chains that don't have a chance of being netted")
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input chain file",
        ))
        .arg(crate::cmd_pgr::args::chain_t_sizes_arg())
        .arg(crate::cmd_pgr::args::chain_q_sizes_arg())
        .arg(crate::cmd_pgr::args::outfile_arg_required())
        .arg(
            Arg::new("dots")
                .long("dots")
                .value_parser(clap::value_parser!(usize))
                .help("Output a dot every so often"),
        )
        .arg(
            Arg::new("pad")
                .long("pad")
                .default_value("1")
                .value_parser(clap::value_parser!(u64))
                .help("Extra to pad around blocks to decrease trash"),
        )
        .arg(crate::cmd_pgr::args::incl_hap_arg())
}
/// Execute the pre-net command.
pub fn execute(args: &ArgMatches) -> Result<()> {
    let input_path = args.get_one::<String>("infile").unwrap();
    let target_sizes_path = args.get_one::<String>("t_sizes").unwrap();
    let query_sizes_path = args.get_one::<String>("q_sizes").unwrap();
    let output_path = crate::cmd_pgr::args::get_outfile(args);

    let dots = args.get_one::<usize>("dots").copied();
    let pad = args.get_one::<u64>("pad").copied().unwrap_or(1);
    let incl_hap = args.get_flag("incl_hap");

    let mut t_hash: HashMap<String, BitMap> = pgr::read_sizes::<u64>(target_sizes_path)?
        .into_iter()
        .map(|(k, v)| (k, BitMap::new(v)))
        .collect();
    let mut q_hash: HashMap<String, BitMap> = pgr::read_sizes::<u64>(query_sizes_path)?
        .into_iter()
        .map(|(k, v)| (k, BitMap::new(v)))
        .collect();

    let reader = pgr::reader(input_path)?;
    let writer = pgr::writer(output_path)?;
    let opts = pgr::libs::chain::PreNetOptions {
        pad,
        incl_hap,
        dots,
    };
    pgr::libs::chain::pre_net(reader, writer, &mut t_hash, &mut q_hash, &opts)
}
