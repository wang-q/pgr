use anyhow::Context;
use clap::{Arg, ArgMatches, Command};

use pgr::libs::fmt::psl::{PslStatsMode, PslStatsOptions};
use std::io::Write;
/// Build the clap subcommand for stats.
pub fn make_subcommand() -> Command {
    Command::new("stats")
        .about("Collects statistics from a PSL file")
        .after_help(
            r###"
Collect statistics from a PSL file.

Examples:
1. Output per-alignment statistics:
   pgr psl stats in.psl -o out.stats

2. Output per-query statistics:
   pgr psl stats --query-stats in.psl -o out.stats

3. Output overall statistics:
   pgr psl stats --overall-stats in.psl -o out.stats
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg().help("Input PSL file. [stdin] for standard input"))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("query_stats")
                .long("query-stats")
                .action(clap::ArgAction::SetTrue)
                .help("Output per-query statistics, the default is per-alignment stats")
                .conflicts_with("overall_stats"),
        )
        .arg(
            Arg::new("overall_stats")
                .long("overall-stats")
                .action(clap::ArgAction::SetTrue)
                .help("Output overall statistics")
                .conflicts_with("query_stats"),
        )
        .arg(
            Arg::new("queries")
                .long("queries")
                .help("Tab separated file with expected qNames and sizes. If specified, statistic will include queries that didn't align."),
        )
        .arg(
            Arg::new("tsv")
                .long("tsv")
                .action(clap::ArgAction::SetTrue)
                .help("Write a TSV header instead of an autoSql header"),
        )
}
/// Execute the stats command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = crate::cmd_pgr::args::get_infile(args);
    let output = crate::cmd_pgr::args::get_outfile(args);
    let query_stats = args.get_flag("query_stats");
    let overall_stats = args.get_flag("overall_stats");
    let queries_file = args.get_one::<String>("queries");
    let tsv = args.get_flag("tsv");

    let reader =
        pgr::reader(input).with_context(|| format!("Failed to open reader for {}", input))?;
    let mut writer =
        pgr::writer(output).with_context(|| format!("Failed to open writer for {}", output))?;

    let mode = if query_stats {
        PslStatsMode::PerQuery
    } else if overall_stats {
        PslStatsMode::Overall
    } else {
        PslStatsMode::PerAlignment
    };

    let opts = PslStatsOptions { mode, tsv };

    let queries = if let Some(q_file) = queries_file {
        let q_reader =
            pgr::reader(q_file).with_context(|| format!("Failed to open reader for {}", q_file))?;
        Some(pgr::libs::fmt::psl::read_queries(q_reader)?)
    } else {
        None
    };

    pgr::libs::fmt::psl::run_stats(reader, &mut writer, &opts, queries)?;

    writer.flush()?;
    Ok(())
}
