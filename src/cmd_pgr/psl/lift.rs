use clap::*;
use std::collections::HashMap;
use std::io::{BufRead, Write};
use std::str::FromStr;

use pgr::libs::psl::Psl;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("lift")
        .about("Lift PSL coordinates from query fragments (e.g., chr1:100-200) to genomic coordinates")
        .after_help(
            r###"
Lift PSL coordinates from query fragments to genomic coordinates.

Notes:
* The query or target name must be in the format `chr:start-end`.
* The coordinates in the name are 1-based, inclusive (UCSC format).
* Requires a chromosome sizes file for correct negative strand lifting.

Examples:
1. Lift query coordinates:
   pgr psl lift input.psl --q-sizes chrom.sizes > output.psl

2. Lift target coordinates:
   pgr psl lift input.psl --t-sizes chrom.sizes > output.psl

3. Lift both:
   pgr psl lift input.psl --q-sizes q.sizes --t-sizes t.sizes > output.psl
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
            Arg::new("q_sizes")
                .long("q-sizes")
                .num_args(1)
                .help("File containing query sizes (name <tab> size)"),
        )
        .arg(
            Arg::new("t_sizes")
                .long("t-sizes")
                .num_args(1)
                .help("File containing target sizes (name <tab> size)"),
        )
}

fn parse_subrange(name: &str) -> Option<(String, u32, u32)> {
    let rg = intspan::Range::from_str(name);
    if rg.is_valid() {
        return Some((rg.chr().to_string(), *rg.start() as u32, *rg.end() as u32));
    }
    None
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap());
    let infile = args.get_one::<String>("infile").unwrap();
    let reader = pgr::reader(infile);

    let q_sizes_file = args.get_one::<String>("q_sizes").map(|s| s.as_str());
    let t_sizes_file = args.get_one::<String>("t_sizes").map(|s| s.as_str());

    // Load sizes if provided
    let load_sizes = |path: Option<&str>| -> HashMap<String, u32> {
        let mut map = HashMap::new();
        if let Some(p) = path {
            let reader = pgr::reader(p);
            for line in reader.lines() {
                if let Ok(line) = line {
                    let parts: Vec<&str> = line.split('\t').collect();
                    if parts.len() >= 2 {
                        if let Ok(size) = parts[1].parse::<u32>() {
                            map.insert(parts[0].to_string(), size);
                        }
                    }
                }
            }
        }
        map
    };

    let q_sizes_map = if q_sizes_file.is_some() { Some(load_sizes(q_sizes_file)) } else { None };
    let t_sizes_map = if t_sizes_file.is_some() { Some(load_sizes(t_sizes_file)) } else { None };

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            if Psl::from_str(&line).is_err() {
                writer.write_fmt(format_args!("{}\n", line))?;
                continue;
            }
        }

        let mut psl: Psl = match line.parse() {
            Ok(p) => p,
            Err(_) => {
                writer.write_fmt(format_args!("{}\n", line))?;
                continue;
            }
        };

        // Try to lift query
        if let Some(sizes_map) = &q_sizes_map {
            if let Some((name_part, start, end)) = parse_subrange(&psl.q_name) {
                let start_0based = start.saturating_sub(1);
                let end_0based = end;
                
                // Check sizes
                let real_size_opt = sizes_map.get(&name_part).copied();
                
                // If sizes provided and match, proceed
                if let Some(real_size) = real_size_opt {
                    if end_0based > real_size {
                        eprintln!("Warning: Subrange end {} > sequence size {} for {}. Skipping query lift.", end_0based, real_size, psl.q_name);
                    } else {
                        let is_neg = psl.strand.starts_with('-');
                        
                        psl.q_name = name_part;
                        psl.q_size = real_size;
                        
                        let offset = if is_neg {
                            real_size - end_0based
                        } else {
                            start_0based
                        };

                        psl.q_start = (psl.q_start as u32 + offset) as i32;
                        psl.q_end = (psl.q_end as u32 + offset) as i32;
                        for q_start in &mut psl.q_starts {
                            *q_start += offset;
                        }
                    }
                } else {
                     eprintln!("Warning: No sizes provided for {}. Skipping query lift.", name_part);
                }
            }
        }

        // Try to lift target
        if let Some(sizes_map) = &t_sizes_map {
            if let Some((name_part, start, end)) = parse_subrange(&psl.t_name) {
                 let start_0based = start.saturating_sub(1);
                let end_0based = end;
                
                // Check sizes
                let real_size_opt = sizes_map.get(&name_part).copied();
                
                // If sizes provided and match, proceed
                if let Some(real_size) = real_size_opt {
                    if end_0based > real_size {
                        eprintln!("Warning: Subrange end {} > sequence size {} for {}. Skipping target lift.", end_0based, real_size, psl.t_name);
                    } else {
                        // For target, check strand if present
                        let is_neg = if psl.strand.len() >= 2 {
                            psl.strand.chars().nth(1).unwrap() == '-'
                        } else {
                            false
                        };

                        psl.t_name = name_part;
                        psl.t_size = real_size;
                        
                        let offset = if is_neg {
                            real_size - end_0based
                        } else {
                            start_0based
                        };

                        psl.t_start = (psl.t_start as u32 + offset) as i32;
                        psl.t_end = (psl.t_end as u32 + offset) as i32;
                        for t_start in &mut psl.t_starts {
                            *t_start += offset;
                        }
                    }
                } else {
                     eprintln!("Warning: No sizes provided for {}. Skipping target lift.", name_part);
                }
            }
        }

        psl.write_to(&mut writer)?;
    }

    Ok(())
}
