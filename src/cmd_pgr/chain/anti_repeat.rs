use anyhow::Context;
use clap::{Arg, ArgMatches, Command};

use pgr::libs::chain::anti_repeat::check_chain;
use pgr::libs::chain::ChainReader;
use pgr::libs::fmt::twobit::TwoBitFile;
use std::io::Write;
// Default scores from UCSC chainAntiRepeat.c
/// Build the clap subcommand for anti-repeat.
pub fn make_subcommand() -> Command {
    Command::new("anti-repeat")
        .about("Filters chains for repeats and degeneracy")
        .after_help(
            r###"
Filters chains for repeats and degeneracy, mirroring the UCSC chainAntiRepeat
workflow. Chains below `--min-score` are dropped; chains above
`--no-check-score` are passed through unchecked; intermediate chains are
evaluated against the target/query 2bit sequence libraries.

Notes:
* Target and query inputs must be 2bit files (use `pgr 2bit` to convert)
* `--min-score` (default: 5000) drops chains below this score
* `--no-check-score` (default: 200000) skips checks above this score

Examples:
1. Basic repeat filtering:
   pgr chain anti-repeat --target-2bit t.2bit --query-2bit q.2bit in.chain -o out.chain

2. Raise the minimum score threshold:
   pgr chain anti-repeat --target-2bit t.2bit --query-2bit q.2bit in.chain --min-score 10000 -o out.chain

3. Lower the no-check threshold:
   pgr chain anti-repeat --target-2bit t.2bit --query-2bit q.2bit in.chain --no-check-score 100000 -o out.chain

"###,
        )
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
/// Execute the anti-repeat command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let target_path = args.get_one::<String>("target_2bit").unwrap();
    let query_path = args.get_one::<String>("query_2bit").unwrap();
    let input_path = args.get_one::<String>("infile").unwrap();
    let output_path = crate::cmd_pgr::args::get_outfile(args);
    let min_score = *args.get_one::<f64>("min_score").unwrap();
    let no_check_score = *args.get_one::<usize>("no_check_score").unwrap();

    let mut target_2bit = TwoBitFile::open(target_path)
        .with_context(|| format!("Failed to open 2bit file {}", target_path))?;
    let mut query_2bit = TwoBitFile::open(query_path)
        .with_context(|| format!("Failed to open 2bit file {}", query_path))?;

    let mut writer = pgr::writer(output_path)
        .with_context(|| format!("Failed to open writer for {}", output_path))?;

    let mut reader = ChainReader::new(
        pgr::reader(input_path)
            .with_context(|| format!("Failed to open reader for {}", input_path))?,
    );

    for res in reader.by_ref() {
        let chain = res?;
        if chain.header.score >= no_check_score as f64 {
            chain.write(&mut writer)?;
            continue;
        }

        if check_chain(&chain, &mut target_2bit, &mut query_2bit, min_score) {
            chain.write(&mut writer)?;
        }
    }

    writer.flush()?;
    Ok(())
}
