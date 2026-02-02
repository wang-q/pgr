use clap::{Arg, ArgMatches, Command};
use std::collections::HashMap;
use std::io::{BufRead, Write};

use pgr::libs::psl::Psl;

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
        .arg(
                    Arg::new("input")
                        .help("Input PSL file")
                        .default_value("stdin")
                        .index(1),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Output file")
                        .default_value("stdout"),
                )
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

#[derive(Debug, Clone, Default)]
struct SumStats {
    q_name: String,
    query_cnt: u32,
    min_q_size: u32,
    max_q_size: u32,
    total_q_size: u64,
    aln_cnt: u32,
    total_align: u64,
    total_match: u64,
    total_rep_match: u64,
    min_ident: f32,
    max_ident: f32,
    min_q_cover: f32,
    max_q_cover: f32,
    min_t_cover: f32,
    max_t_cover: f32,
    min_rep_match: f32,
    max_rep_match: f32,
}

impl SumStats {
    fn new(q_name: &str, q_size: u32) -> Self {
        Self {
            q_name: q_name.to_string(),
            query_cnt: 1,
            min_q_size: q_size,
            max_q_size: q_size,
            total_q_size: 0,
            aln_cnt: 0,
            ..Default::default()
        }
    }

    fn accumulate(&mut self, psl: &Psl) {
        let ident = psl.calc_ident();
        let q_cover = psl.calc_q_cover();
        let t_cover = psl.calc_t_cover();
        let rep_match = psl.calc_rep_match();

        self.total_q_size += psl.q_size as u64;

        if self.aln_cnt == 0 {
            self.min_ident = ident;
            self.max_ident = ident;
            self.min_q_cover = q_cover;
            self.max_q_cover = q_cover;
            self.min_t_cover = t_cover;
            self.max_t_cover = t_cover;
            self.min_rep_match = rep_match;
            self.max_rep_match = rep_match;

            self.min_q_size = self.min_q_size.min(psl.q_size);
            self.max_q_size = self.max_q_size.max(psl.q_size);
        } else {
            self.min_q_size = self.min_q_size.min(psl.q_size);
            self.max_q_size = self.max_q_size.max(psl.q_size);

            self.min_ident = self.min_ident.min(ident);
            self.max_ident = self.max_ident.max(ident);

            self.min_q_cover = self.min_q_cover.min(q_cover);
            self.max_q_cover = self.max_q_cover.max(q_cover);

            self.min_t_cover = self.min_t_cover.min(t_cover);
            self.max_t_cover = self.max_t_cover.max(t_cover);

            self.min_rep_match = self.min_rep_match.min(rep_match);
            self.max_rep_match = self.max_rep_match.max(rep_match);
        }

        self.total_align += psl.calc_aligned() as u64;
        self.total_match += psl.calc_match() as u64;
        self.total_rep_match += psl.rep_match as u64;
        self.aln_cnt += 1;
    }

    fn merge(&mut self, other: &SumStats) {
        if self.aln_cnt == 0 {
            self.min_q_size = other.min_q_size;
            self.max_q_size = other.max_q_size;
            self.min_ident = other.min_ident;
            self.max_ident = other.max_ident;
            self.min_q_cover = other.min_q_cover;
            self.max_q_cover = other.max_q_cover;
            self.min_t_cover = other.min_t_cover;
            self.max_t_cover = other.max_t_cover;
            self.min_rep_match = other.min_rep_match;
            self.max_rep_match = other.max_rep_match;
        } else if other.aln_cnt > 0 {
            self.min_q_size = self.min_q_size.min(other.min_q_size);
            self.max_q_size = self.max_q_size.max(other.max_q_size);
            self.min_ident = self.min_ident.min(other.min_ident);
            self.max_ident = self.max_ident.max(other.max_ident);
            self.min_q_cover = self.min_q_cover.min(other.min_q_cover);
            self.max_q_cover = self.max_q_cover.max(other.max_q_cover);
            self.min_t_cover = self.min_t_cover.min(other.min_t_cover);
            self.max_t_cover = self.max_t_cover.max(other.max_t_cover);
            self.min_rep_match = self.min_rep_match.min(other.min_rep_match);
            self.max_rep_match = self.max_rep_match.max(other.max_rep_match);
        }

        self.query_cnt += other.query_cnt;
        self.total_q_size += other.total_q_size;
        self.total_align += other.total_align;
        self.total_match += other.total_match;
        self.total_rep_match += other.total_rep_match;
        self.aln_cnt += other.aln_cnt;
    }

    fn mean_ident(&self) -> f32 {
        if self.total_align == 0 {
            0.0
        } else {
            self.total_match as f32 / self.total_align as f32
        }
    }

    fn mean_q_size(&self) -> u32 {
        if self.query_cnt == 0 {
            0
        } else {
            (self.total_q_size / self.query_cnt as u64) as u32
        }
    }

    fn mean_q_cover(&self) -> f32 {
        if self.total_q_size == 0 {
            0.0
        } else {
            self.total_align as f32 / self.total_q_size as f32
        }
    }

    fn mean_rep_match(&self) -> f32 {
        if self.total_align == 0 {
            0.0
        } else {
            self.total_rep_match as f32 / self.total_align as f32
        }
    }
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = args.get_one::<String>("input").unwrap();
    let output = args.get_one::<String>("output").unwrap();
    let query_stats = args.get_flag("query_stats");
    let overall_stats = args.get_flag("overall_stats");
    let queries_file = args.get_one::<String>("queries");
    let tsv = args.get_flag("tsv");

    let reader = intspan::reader(input);
    let mut writer = intspan::writer(output);

    let mut query_stats_tbl: HashMap<String, SumStats> = HashMap::new();

    // Load queries if provided
    if let Some(q_file) = queries_file {
        let r = intspan::reader(q_file);
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
        for line in reader.lines() {
            let line = line?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            if let Ok(psl) = line.parse::<Psl>() {
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

        for line in reader.lines() {
            let line = line?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            if let Ok(psl) = line.parse::<Psl>() {
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
