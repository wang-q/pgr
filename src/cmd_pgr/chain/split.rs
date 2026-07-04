use clap::{Arg, ArgMatches, Command};
use pgr::libs::chain::ChainReader;
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::Path;
/// Build the clap subcommand for split.
pub fn make_subcommand() -> Command {
    Command::new("split")
        .about("Splits chains up by target or query sequence")
        .arg(crate::cmd_pgr::args::infiles_arg("chain"))
        .arg(crate::cmd_pgr::args::outdir_arg_required())
        .arg(crate::cmd_pgr::args::by_query_arg(
            "Split on query (default is on target)",
        ))
        .arg(
            Arg::new("lump")
                .long("lump")
                .value_parser(clap::value_parser!(usize))
                .help("Lump together so have only N split files"),
        )
}
/// Execute the split command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let out_dir = args.get_one::<String>("outdir").unwrap();
    let chain_files: Vec<_> = args.get_many::<String>("infiles").unwrap().collect();
    let split_on_q = args.get_flag("by_query");
    let lump = args.get_one::<usize>("lump").copied();

    fs::create_dir_all(out_dir)?;

    // Cache open file handles
    let mut file_cache: HashMap<String, Box<dyn Write>> = HashMap::new();

    for file_path in chain_files {
        let reader = ChainReader::new(pgr::reader(file_path)?);

        for res in reader {
            let chain = res?;

            let raw_name = if split_on_q {
                &chain.header.q_name
            } else {
                &chain.header.t_name
            };

            let name = if let Some(l) = lump {
                pgr::libs::chain::lump_name(raw_name, l)
            } else {
                raw_name.clone()
            };

            let writer = match file_cache.entry(name.clone()) {
                std::collections::hash_map::Entry::Occupied(e) => e.into_mut(),
                std::collections::hash_map::Entry::Vacant(e) => {
                    let path = Path::new(out_dir).join(format!("{}.chain", name));
                    // Truncate on first open in this process; subsequent writes via cache.
                    let path_str = path
                        .to_str()
                        .ok_or_else(|| anyhow::anyhow!("non-UTF-8 path"))?;
                    e.insert(pgr::writer(path_str)?)
                }
            };

            chain.write(writer)?;
        }
    }

    Ok(())
}
