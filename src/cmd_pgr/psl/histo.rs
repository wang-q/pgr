use clap::*;
use pgr::libs::fmt::psl::Psl;
use std::collections::HashMap;

pub fn make_subcommand() -> Command {
    Command::new("histo")
        .about("Collect counts on PSL alignments for making histograms")
        .after_help(
            r###"
These then be analyzed with R, textHistogram, etc.

The 'what' argument determines what data to collect, the following are currently supported:

* alignsPerQuery - number of alignments per query. Output is one line per query with the number of alignments.

* coverSpread - difference between the highest and lowest coverage for alignments of a query. Output line per query, with the difference. Use --multi-only to omit queries with a single alignment.

* idSpread - difference between the highest and lowest fraction identity for alignments of a query. Output line per query, with the difference.

Examples:
  # Collect alignment counts per query
  pgr psl histo --what alignsPerQuery in.psl -o out.histo
"###,
        )
        .arg(
            Arg::new("what")
                .long("what")
                .required(true)
                .value_name("TYPE")
                .value_parser(["alignsPerQuery", "coverSpread", "idSpread"])
                .help("What data to collect"),
        )
        .arg(crate::cmd_pgr::args::infile_arg().help("Input PSL file"))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("multi_only")
                .long("multi-only")
                .short('m')
                .action(ArgAction::SetTrue)
                .help("Omit queries with only one alignment"),
        )
        .arg(
            Arg::new("non_zero")
                .long("non-zero")
                .short('z')
                .action(ArgAction::SetTrue)
                .help("Omit queries with zero values"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let what = args.get_one::<String>("what").unwrap();
    let input = crate::cmd_pgr::args::get_infile(args);
    let output = crate::cmd_pgr::args::get_outfile(args);
    let multi_only = args.get_flag("multi_only");
    let non_zero = args.get_flag("non_zero");

    let reader = pgr::reader(input)?;
    let mut writer = pgr::writer(output)?;

    // Read all PSLs and group by query
    let mut query_map: HashMap<String, Vec<Psl>> = HashMap::new();

    for psl in pgr::libs::fmt::psl::iter_psl(reader) {
        let psl = psl?;
        query_map.entry(psl.q_name.clone()).or_default().push(psl);
    }

    // Process queries (iteration order not guaranteed, but usually fine for histograms.
    // If output order matters, we should sort keys. C implementation uses hash table, likely random order.)
    // Let's sort keys for deterministic output.
    let mut queries: Vec<_> = query_map.keys().cloned().collect();
    queries.sort();

    for q_name in queries {
        let psls = &query_map[&q_name];

        if multi_only && psls.len() <= 1 {
            continue;
        }

        match what.as_str() {
            "alignsPerQuery" => {
                let cnt = psls.len();
                if !non_zero || cnt != 0 {
                    // cnt is never 0 here if it exists in map, but logic follows C
                    writeln!(writer, "{}", cnt)?;
                }
            }
            "coverSpread" => {
                let (min, max) = pgr::libs::fmt::psl::calc_spread(psls, |p| p.cover());
                let diff = max - min;
                if !non_zero || diff != 0.0 {
                    // Using same format as C: %.4g
                    // Rust doesn't have %g exactly, but {:.*} might work or standard formatting.
                    // C uses %0.4g.
                    // Let's use generic formatting for now, maybe check precision requirements.
                    // %g uses scientific notation for large/small numbers.
                    writeln!(writer, "{:.4}", diff)?;
                }
            }
            "idSpread" => {
                let (min, max) = pgr::libs::fmt::psl::calc_spread(psls, |p| p.ident());
                let diff = max - min;
                if !non_zero || diff != 0.0 {
                    writeln!(writer, "{:.4}", diff)?;
                }
            }
            _ => anyhow::bail!("unsupported stat type"),
        }
    }

    Ok(())
}
