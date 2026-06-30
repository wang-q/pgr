use clap::{Arg, ArgMatches, Command};
use std::collections::HashMap;
use std::io::{BufRead, Write};

use pgr::libs::fmt::psl::SumStats;

pub fn make_subcommand() -> Command {
    Command::new("stats")
        .about("Collect statistics from a psl file")
        .after_help(
            r###"
Collect statistics from a psl file.

Examples:
  pgr psl stats in.psl -o out.stats
  pgr psl stats --query-stats in.psl -o out.stats
  pgr psl stats --overall-stats in.psl -o out.stats
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg().help("Input PSL file"))
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

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = crate::cmd_pgr::args::get_infile(args);
    let output = crate::cmd_pgr::args::get_outfile(args);
    let query_stats = args.get_flag("query_stats");
    let overall_stats = args.get_flag("overall_stats");
    let queries_file = args.get_one::<String>("queries");
    let tsv = args.get_flag("tsv");

    let reader = pgr::reader(input)?;
    let mut writer = pgr::writer(output)?;

    let mut query_stats_tbl: HashMap<String, SumStats> = HashMap::new();

    // Load queries if provided
    if let Some(q_file) = queries_file {
        let r = pgr::reader(q_file)?;
        for line in r.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 2 {
                let q_name = parts[0].to_string();
                let q_size: u32 = parts[1].parse()?;
                query_stats_tbl.insert(q_name.clone(), SumStats::new(&q_name, q_size));
            }
        }
    }

    if query_stats || overall_stats {
        // Aggregation modes
        for psl in pgr::libs::fmt::psl::iter_psl(reader) {
            let psl = psl?;
            if queries_file.is_some() {
                if let Some(entry) = query_stats_tbl.get_mut(&psl.q_name) {
                    entry.accumulate(&psl);
                }
            } else {
                let entry = query_stats_tbl
                    .entry(psl.q_name.clone())
                    .or_insert_with(|| SumStats::new(&psl.q_name, psl.q_size));
                entry.accumulate(&psl);
            }
        }

        if query_stats {
            if !tsv {
                write!(writer, "#")?;
            }
            writeln!(writer, "qName\tqSize\talnCnt\tminIdent\tmaxIdent\tmeanIdent\tminQCover\tmaxQCover\tmeanQCover\tminRepMatch\tmaxRepMatch\tmeanRepMatch\tminTCover\tmaxTCover")?;

            let mut keys: Vec<_> = query_stats_tbl.keys().cloned().collect();
            keys.sort();

            for k in keys {
                let s = &query_stats_tbl[&k];
                writeln!(writer, "{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
                    s.q_name, s.min_q_size, s.aln_cnt,
                    s.min_ident, s.max_ident, s.mean_ident(),
                    s.min_q_cover, s.max_q_cover, s.mean_q_cover(),
                    s.min_rep_match, s.max_rep_match, s.mean_rep_match(),
                    s.min_t_cover, s.max_t_cover
                )?;
            }
        } else {
            // overall stats
            let mut os = SumStats::default();
            let mut aligned1 = 0;
            let mut aligned_n = 0;

            for s in query_stats_tbl.values() {
                os.merge(s);

                if s.aln_cnt == 1 {
                    aligned1 += 1;
                } else if s.aln_cnt > 1 {
                    aligned_n += 1;
                }
            }

            if !tsv {
                write!(writer, "#")?;
            }
            writeln!(writer, "queryCnt\tminQSize\tmaxQSize\tmeanQSize\talnCnt\tminIdent\tmaxIdent\tmeanIdent\tminQCover\tmaxQCover\tmeanQCover\tminRepMatch\tmaxRepMatch\tmeanRepMatch\tminTCover\tmaxTCover\taligned\taligned1\talignedN\ttotalAlignedSize")?;

            writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}",
                os.query_cnt, os.min_q_size, os.max_q_size, os.mean_q_size(),
                os.aln_cnt,
                os.min_ident, os.max_ident, os.mean_ident(),
                os.min_q_cover, os.max_q_cover, os.mean_q_cover(),
                os.min_rep_match, os.max_rep_match, os.mean_rep_match(),
                os.min_t_cover, os.max_t_cover,
                aligned1 + aligned_n, aligned1, aligned_n,
                os.total_align
            )?;
        }
    } else {
        // Per-alignment mode
        if !tsv {
            write!(writer, "#")?;
        }
        writeln!(
            writer,
            "qName\tqSize\ttName\ttStart\ttEnd\tident\tqCover\trepMatch\ttCover"
        )?;

        for psl in pgr::libs::fmt::psl::iter_psl(reader) {
            let psl = psl?;
            if queries_file.is_some() {
                if let Some(entry) = query_stats_tbl.get_mut(&psl.q_name) {
                    writeln!(
                        writer,
                        "{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
                        psl.q_name,
                        psl.q_size,
                        psl.t_name,
                        psl.t_start,
                        psl.t_end,
                        psl.calc_ident(),
                        psl.calc_q_cover(),
                        psl.calc_rep_match(),
                        psl.calc_t_cover()
                    )?;
                    entry.aln_cnt += 1;
                }
            } else {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
                    psl.q_name,
                    psl.q_size,
                    psl.t_name,
                    psl.t_start,
                    psl.t_end,
                    psl.calc_ident(),
                    psl.calc_q_cover(),
                    psl.calc_rep_match(),
                    psl.calc_t_cover()
                )?;
            }
        }

        if queries_file.is_some() {
            let mut keys: Vec<_> = query_stats_tbl.keys().cloned().collect();
            keys.sort();

            for k in keys {
                let s = &query_stats_tbl[&k];
                if s.aln_cnt == 0 {
                    writeln!(
                        writer,
                        "{}\t{}\t\t0\t0\t0.0000\t0.0000\t0.0000\t0.0000",
                        s.q_name, s.min_q_size
                    )?;
                }
            }
        }
    }

    Ok(())
}
