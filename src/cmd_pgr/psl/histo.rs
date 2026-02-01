use clap::*;
use pgr::libs::psl::Psl;
use std::collections::HashMap;
use std::io::BufRead;
use std::str::FromStr;

pub fn make_subcommand() -> Command {
    Command::new("histo")
        .about("Collect counts on PSL alignments for making histograms")
        .arg(
            Arg::new("what")
                .long("what")
                .required(true)
                .value_parser(["alignsPerQuery", "coverSpread", "idSpread"])
                .help("What data to collect"),
        )
        .arg(
            Arg::new("input")
                .index(1)
                .help("Input PSL file")
                .default_value("stdin"),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .help("Output histogram file")
                .default_value("stdout"),
        )
        .arg(
            Arg::new("multi_only")
                .long("multi-only")
                .action(ArgAction::SetTrue)
                .help("Omit queries with only one alignment"),
        )
        .arg(
            Arg::new("non_zero")
                .long("non-zero")
                .action(ArgAction::SetTrue)
                .help("Omit queries with zero values"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let what = args.get_one::<String>("what").unwrap();
    let input = args.get_one::<String>("input").unwrap();
    let output = args.get_one::<String>("output").unwrap();
    let multi_only = args.get_flag("multi_only");
    let non_zero = args.get_flag("non_zero");

    let reader = intspan::reader(input);
    let mut writer = intspan::writer(output);

    // Read all PSLs and group by query
    let mut query_map: HashMap<String, Vec<Psl>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        
        // Skip header lines often found in PSL files
        // psLayout version 3
        // match ...
        // ---------------- ...
        if line.starts_with("psLayout") || line.starts_with("match") || line.starts_with("------") {
             continue;
        }

        let psl = Psl::from_str(&line)?;
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
                let (min, max) = calc_spread(psls, calc_cover);
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
                let (min, max) = calc_spread(psls, calc_ident);
                let diff = max - min;
                if !non_zero || diff != 0.0 {
                    writeln!(writer, "{:.4}", diff)?;
                }
            }
            _ => unreachable!(),
        }
    }

    Ok(())
}

fn calc_cover(psl: &Psl) -> f32 {
    let aligned = psl.match_count + psl.mismatch_count + psl.rep_match;
    if aligned == 0 {
        0.0
    } else {
        aligned as f32 / psl.q_size as f32
    }
}

fn calc_ident(psl: &Psl) -> f32 {
    let aligned = psl.match_count + psl.mismatch_count + psl.rep_match;
    if aligned == 0 {
        0.0
    } else {
        (psl.match_count + psl.rep_match) as f32 / aligned as f32
    }
}

fn calc_spread<F>(psls: &[Psl], func: F) -> (f32, f32)
where
    F: Fn(&Psl) -> f32,
{
    let mut min_val = f32::MAX;
    let mut max_val = f32::MIN;

    for psl in psls {
        let val = func(psl);
        if val < min_val {
            min_val = val;
        }
        if val > max_val {
            max_val = val;
        }
    }
    
    // Handle case where psls is empty (shouldn't happen here)
    if min_val == f32::MAX {
        (0.0, 0.0)
    } else {
        (min_val, max_val)
    }
}
