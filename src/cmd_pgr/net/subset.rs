use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::chain::{read_chains, Chain};
use pgr::libs::net::{read_nets, Fill};
use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::rc::Rc;

pub fn make_subcommand() -> Command {
    Command::new("subset")
        .about("Create chain file with subset of chains that appear in the net")
        .arg(
            Arg::new("net_in")
                .required(true)
                .index(1)
                .help("Input net file"),
        )
        .arg(
            Arg::new("chain_in")
                .required(true)
                .index(2)
                .help("Input chain file"),
        )
        .arg(
            Arg::new("chain_out")
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
        .arg(
            Arg::new("type")
                .long("type")
                .action(ArgAction::Set)
                .help("Restrict output to particular type in net file"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let net_in = args.get_one::<String>("net_in").unwrap();
    let chain_in = args.get_one::<String>("chain_in").unwrap();
    let chain_out = args.get_one::<String>("chain_out").unwrap();
    let whole_chains = args.get_flag("whole_chains");
    let split_on_insert = args.get_flag("split_on_insert");
    let type_filter = args.get_one::<String>("type");

    // Read chains
    let chain_reader = BufReader::new(File::open(chain_in)?);
    let chains_vec = read_chains(chain_reader)?;
    let mut chains_map: HashMap<u64, Chain> = HashMap::new();
    for chain in chains_vec {
        chains_map.insert(chain.header.id, chain);
    }

    // Read nets
    let net_reader = BufReader::new(File::open(net_in)?);
    let chroms = read_nets(net_reader)?;

    let mut writer = BufWriter::new(File::create(chain_out)?);

    for chrom in chroms {
        // Traverse net structure
        // Root is a Gap
        process_gap(
            &chrom.root,
            &chains_map,
            &mut writer,
            whole_chains,
            split_on_insert,
            type_filter,
        )?;
    }

    Ok(())
}

fn process_gap(
    gap: &Rc<RefCell<pgr::libs::net::Gap>>,
    chains_map: &HashMap<u64, Chain>,
    writer: &mut impl std::io::Write,
    whole_chains: bool,
    split_on_insert: bool,
    type_filter: Option<&String>,
) -> anyhow::Result<()> {
    let gap = gap.borrow();
    for fill in &gap.fills {
        process_fill(
            fill,
            chains_map,
            writer,
            whole_chains,
            split_on_insert,
            type_filter,
        )?;
    }
    Ok(())
}

fn process_fill(
    fill_rc: &Rc<RefCell<Fill>>,
    chains_map: &HashMap<u64, Chain>,
    writer: &mut impl std::io::Write,
    whole_chains: bool,
    split_on_insert: bool,
    type_filter: Option<&String>,
) -> anyhow::Result<()> {
    let fill = fill_rc.borrow();

    // Check type filter
    if let Some(t) = type_filter {
        if &fill.class != t {
            return Ok(()); // Skip but continue traversal?
                           // In C: if (!sameString(type, fill->type)) return;
                           // It returns from convertFill, but then it continues recursion in rConvert.
                           // Wait, in C rConvert calls convertFill THEN recurses.
                           // So if type doesn't match, we don't output this fill, but do we recurse?
                           // C code:
                           // if (fill->chainId) { ... convertFill ... }
                           // if (fill->children) rConvert(...);
                           //
                           // convertFill checks type and returns if mismatch.
                           // So yes, we should still recurse.
        }
    }

    // Process current fill
    if fill.chain_id != 0 {
        if let Some(chain) = chains_map.get(&fill.chain_id) {
            if whole_chains {
                chain.write(writer)?;
            } else if split_on_insert {
                // Split on insert logic
                let mut t_start = fill.start;

                // Iterate over gaps to find inserts
                for gap_rc in &fill.gaps {
                    let gap = gap_rc.borrow();
                    if !gap.fills.is_empty() {
                        // This gap has inserts (children fills)
                        // Output chain part from t_start to gap.start
                        if gap.start > t_start {
                            if let Some(sub) = chain.subset(t_start, gap.start) {
                                sub.write(writer)?;
                            }
                        }
                        t_start = gap.end;
                    }
                }
                // Output remaining part
                if fill.end > t_start {
                    if let Some(sub) = chain.subset(t_start, fill.end) {
                        sub.write(writer)?;
                    }
                }
            } else {
                // Default: subset to fill range
                if let Some(sub) = chain.subset(fill.start, fill.end) {
                    sub.write(writer)?;
                }
            }
        }
    }

    // Recurse into children gaps
    for gap in &fill.gaps {
        process_gap(
            gap,
            chains_map,
            writer,
            whole_chains,
            split_on_insert,
            type_filter,
        )?;
    }

    Ok(())
}
