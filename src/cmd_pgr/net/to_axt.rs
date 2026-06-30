use clap::{Arg, ArgMatches, Command};
use pgr::libs::chain::net::{net_to_axt, read_nets};
use pgr::libs::chain::sub_matrix::SubMatrix;
use pgr::libs::chain::{Chain, ChainReader};
use pgr::libs::fmt::twobit::TwoBitFile;
use std::collections::HashMap;
use std::fs::File;

pub fn make_subcommand() -> Command {
    Command::new("to-axt")
        .about("Convert net (and chain) to axt")
        .arg(Arg::new("in_net").required(true).help("Input net file"))
        .arg(Arg::new("in_chain").required(true).help("Input chain file"))
        .arg(Arg::new("target").required(true).help("Target 2bit file"))
        .arg(Arg::new("query").required(true).help("Query 2bit file"))
        .arg(Arg::new("out_axt").required(true).help("Output AXT file"))
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let in_net = matches.get_one::<String>("in_net").unwrap();
    let in_chain = matches.get_one::<String>("in_chain").unwrap();
    let target = matches.get_one::<String>("target").unwrap();
    let query = matches.get_one::<String>("query").unwrap();
    let out_axt = matches.get_one::<String>("out_axt").unwrap();

    let mut t_2bit = TwoBitFile::open(target)?;
    let mut q_2bit = TwoBitFile::open(query)?;

    let mut chains: HashMap<u64, Chain> = HashMap::new();
    let chain_reader = ChainReader::new(File::open(in_chain)?);
    for chain_res in chain_reader {
        let chain = chain_res?;
        chains.insert(chain.header.id, chain);
    }

    let reader = pgr::reader(in_net)?;
    let nets = read_nets(reader)?;

    let matrix = SubMatrix::hoxd55();

    let mut writer = pgr::writer(out_axt)?;
    net_to_axt(
        &nets,
        &chains,
        &mut t_2bit,
        &mut q_2bit,
        &matrix,
        &mut writer,
    )?;

    Ok(())
}
