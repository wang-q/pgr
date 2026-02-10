use clap::*;
use pgr::libs::psl::Psl;
use std::io::{BufRead, Write};
use std::str::FromStr;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("to-range")
        .about("Extract coordinates from PSL as ranges (.rg)")
        .after_help(
            r###"
Extract alignment coordinates from PSL files and output in .rg format (chr:start-end).
This is useful for depth calculation with `spanr coverage`.

Notes:
* Coordinates are converted to 1-based inclusive (intspan/UCSC format).
* Supports strand-aware coordinate conversion (outputs positive strand coordinates).
* Outputs one range per alignment block.

Examples:
1. Extract query ranges:
   pgr psl to-range input.psl > query.rg

2. Extract target ranges:
   pgr psl to-range input.psl --target > target.rg
"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Input PSL file. [stdin] for standard input"),
        )
        .arg(
            Arg::new("outfile")
                .short('o')
                .long("outfile")
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
        .arg(
            Arg::new("target")
                .long("target")
                .short('t')
                .action(ArgAction::SetTrue)
                .help("Extract target coordinates instead of query"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap());
    let infile = args.get_one::<String>("infile").unwrap();
    let reader = pgr::reader(infile);
    let extract_target = args.get_flag("target");

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        // Skip header lines
        if line.starts_with("psLayout") || line.starts_with("match") || line.starts_with("------") {
            continue;
        }

        let psl = match Psl::from_str(&line) {
            Ok(p) => p,
            Err(_) => continue,
        };

        let (name, size, starts, is_neg) = if extract_target {
            let is_neg = if psl.strand.len() >= 2 {
                psl.strand.chars().nth(1).unwrap() == '-'
            } else {
                false
            };
            (&psl.t_name, psl.t_size, &psl.t_starts, is_neg)
        } else {
            let is_neg = psl.strand.starts_with('-');
            (&psl.q_name, psl.q_size, &psl.q_starts, is_neg)
        };

        for (i, &start) in starts.iter().enumerate() {
            let len = psl.block_sizes[i];
            let end = start + len; // 0-based exclusive end

            // Convert to 1-based inclusive range on positive strand
            let (final_start, final_end) = if is_neg {
                // Reverse complement coordinates
                // 0-based: [size - end, size - start)
                // 1-based: size - end + 1, size - start
                (size - end + 1, size - start)
            } else {
                // Positive strand coordinates
                // 0-based: [start, end)
                // 1-based: start + 1, end
                (start + 1, end)
            };

            writer.write_fmt(format_args!("{}:{}-{}\n", name, final_start, final_end))?;
        }
    }

    Ok(())
}
