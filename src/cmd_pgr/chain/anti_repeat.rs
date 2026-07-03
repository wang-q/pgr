use clap::{Arg, ArgMatches, Command};

use pgr::libs::chain::anti_repeat::check_chain;
use pgr::libs::chain::read_chains;
use pgr::libs::fmt::twobit::TwoBitFile;
// Default scores from UCSC chainAntiRepeat.c

pub fn make_subcommand() -> Command {
    Command::new("anti-repeat")
        .about("Filter chains for repeats and degeneracy")
        .arg(crate::cmd_pgr::args::target_2bit_arg())
        .arg(crate::cmd_pgr::args::query_2bit_arg())
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input chain file",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg_required())
        .arg(crate::cmd_pgr::args::min_score_arg("5000"))
        .arg(
            Arg::new("no_check_score")
                .long("no-check-score")
                .default_value("200000")
                .value_parser(clap::value_parser!(usize))
                .help("Score above which no checks are performed"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let target_path = args.get_one::<String>("target_2bit").unwrap();
    let query_path = args.get_one::<String>("query_2bit").unwrap();
    let input_path = args.get_one::<String>("infile").unwrap();
    let output_path = crate::cmd_pgr::args::get_outfile(args);
    let min_score = *args.get_one::<f64>("min_score").unwrap();
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
