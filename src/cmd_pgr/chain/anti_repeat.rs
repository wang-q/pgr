use clap::{Arg, ArgMatches, Command};

use pgr::libs::chain::anti_repeat::check_chain;
use pgr::libs::chain::read_chains;
use pgr::libs::fmt::twobit::TwoBitFile;
// Default scores from UCSC chainAntiRepeat.c

pub fn make_subcommand() -> Command {
    Command::new("anti-repeat")
        .about("Filter chains for repeats and degeneracy")
        .arg(
            Arg::new("target")
                .long("target")
                .short('t')
                .required(true)
                .help("Target genome 2bit file"),
        )
        .arg(
            Arg::new("query")
                .long("query")
                .short('q')
                .required(true)
                .help("Query genome 2bit file"),
        )
        .arg(Arg::new("input").required(true).help("Input chain file"))
        .arg(Arg::new("output").required(true).help("Output chain file"))
        .arg(
            Arg::new("min_score")
                .long("min-score")
                .default_value("5000")
                .value_parser(clap::value_parser!(usize))
                .help("Minimum score to pass"),
        )
        .arg(
            Arg::new("no_check_score")
                .long("no-check-score")
                .default_value("200000")
                .value_parser(clap::value_parser!(usize))
                .help("Score above which no checks are performed"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let target_path = args.get_one::<String>("target").unwrap();
    let query_path = args.get_one::<String>("query").unwrap();
    let input_path = args.get_one::<String>("input").unwrap();
    let output_path = args.get_one::<String>("output").unwrap();
    let min_score = *args.get_one::<usize>("min_score").unwrap();
    let no_check_score = *args.get_one::<usize>("no_check_score").unwrap();

    let mut target_2bit = TwoBitFile::open(target_path)?;
    let mut query_2bit = TwoBitFile::open(query_path)?;

    let chains = read_chains(pgr::reader(input_path)?)?; // Note: read_chains reads all chains into memory

    let mut writer = pgr::writer(output_path)?;

    for chain in chains {
        if chain.header.score >= no_check_score as f64 {
            chain.write(&mut writer)?;
            continue;
        }

        if check_chain(&chain, &mut target_2bit, &mut query_2bit, min_score) {
            chain.write(&mut writer)?;
        }
    }

    Ok(())
}
