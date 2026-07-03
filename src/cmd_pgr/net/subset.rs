use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::chain::net::{read_nets, subset_nets, SubsetOptions};
use pgr::libs::chain::{read_chains, Chain};
use std::collections::HashMap;

pub fn make_subcommand() -> Command {
    Command::new("subset")
        .about("Create chain file with subset of chains that appear in the net")
        .arg(crate::cmd_pgr::args::in_net_arg().index(1))
        .arg(crate::cmd_pgr::args::in_chain_arg().index(2))
        .arg(
            Arg::new("out_chain")
                .required(true)
                .index(3)
                .help("Output chain file"),
        )
        .arg(
            Arg::new("whole_chains")
                .long("whole-chains")
                .action(ArgAction::SetTrue)
                .help("Write entire chain references by net, don't split"),
        )
        .arg(
            Arg::new("split_on_insert")
                .long("split-on-insert")
                .action(ArgAction::SetTrue)
                .help("Split chain when get an insertion of another chain"),
        )
        .arg(crate::cmd_pgr::args::net_type_arg(
            ArgAction::Set,
            "Restrict output to particular type in net file",
        ))
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let net_in = args.get_one::<String>("in_net").unwrap();
    let chain_in = args.get_one::<String>("in_chain").unwrap();
    let chain_out = args.get_one::<String>("out_chain").unwrap();
    let whole_chains = args.get_flag("whole_chains");
    let split_on_insert = args.get_flag("split_on_insert");
    let type_filter = args.get_one::<String>("type");

    // Read chains
    let chain_reader = pgr::reader(chain_in)?;
    let chains_vec = read_chains(chain_reader)?;
    let mut chains_map: HashMap<u64, Chain> = HashMap::new();
    for chain in chains_vec {
        chains_map.insert(chain.header.id, chain);
    }

    // Read nets
    let net_reader = pgr::reader(net_in)?;
    let chroms = read_nets(net_reader)?;

    let mut writer = pgr::writer(chain_out)?;

    let opts = SubsetOptions {
        whole_chains,
        split_on_insert,
    };
    subset_nets(&chroms, &chains_map, &mut writer, opts, type_filter)?;

    Ok(())
}
