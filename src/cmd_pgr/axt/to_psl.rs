use clap::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use pgr::libs::axt::AxtReader;
use pgr::libs::psl::Psl;

pub fn make_subcommand() -> Command {
    Command::new("to-psl")
        .about("Convert from axt to psl format")
        .after_help(
            r###"
Where tSizes and qSizes are tab-delimited files with <seqName> <size> columns.

Examples:
  # Convert axt to psl
  pgr axt to-psl in.axt -t t.sizes -q q.sizes -o out.psl
"###,
        )
        .arg(
            Arg::new("input")
                .help("Input AXT file")
                .default_value("stdin")
                .index(1),
        )
        .arg(
            Arg::new("t_sizes")
                .long("t-sizes")
                .short('t')
                .value_name("FILE")
                .help("Target sizes file")
                .required(true),
        )
        .arg(
            Arg::new("q_sizes")
                .long("q-sizes")
                .short('q')
                .value_name("FILE")
                .help("Query sizes file")
                .required(true),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("FILE")
                .help("Output PSL file")
                .default_value("stdout"),
        )
}

fn load_sizes(path: &str) -> anyhow::Result<HashMap<String, usize>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut sizes = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let name = parts[0].to_string();
            let size = parts[1].parse::<usize>()?;
            sizes.insert(name, size);
        }
    }

    Ok(sizes)
}

fn reverse_range(start: &mut i32, end: &mut i32, size: u32) {
    let s = *start;
    let e = *end;
    *start = size as i32 - e;
    *end = size as i32 - s;
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = args.get_one::<String>("input").unwrap();
    let t_sizes_path = args.get_one::<String>("t_sizes").unwrap();
    let q_sizes_path = args.get_one::<String>("q_sizes").unwrap();
    let output = args.get_one::<String>("output").unwrap();

    // Load sizes
    let t_sizes = load_sizes(t_sizes_path)?;
    let q_sizes = load_sizes(q_sizes_path)?;

    // Open input
    let reader = intspan::reader(input);
    let reader = AxtReader::new(reader);

    // Open output
    let mut writer = intspan::writer(output);

    for result in reader {
        let axt = result?;

        // Get sizes
        let q_size = *q_sizes
            .get(&axt.q_name)
            .ok_or_else(|| anyhow::anyhow!("Query size not found for {}", axt.q_name))?;
        let t_size = *t_sizes
            .get(&axt.t_name)
            .ok_or_else(|| anyhow::anyhow!("Target size not found for {}", axt.t_name))?;

        // Prepare coordinates
        // libs/axt.rs returns 0-based half-open coordinates
        let mut q_start = axt.q_start as i32;
        let mut q_end = axt.q_end as i32;

        // axtToPsl.c logic: "if (axt->qStrand == '-') reverseIntRange(&qStart, &qEnd, qSize);"
        // This converts strand-relative coordinates (as in AXT) to positive strand coordinates
        // which pslFromAlign expects (so it can reverse them back internally).
        if axt.q_strand == '-' {
            reverse_range(&mut q_start, &mut q_end, q_size as u32);
        }

        // Construct strand string for PSL (e.g. "-")
        // Note: PSL usually tracks target strand implicitly as +, so strand field is just query strand?
        // axtToPsl.c: strand[0] = axt->qStrand; strand[1] = '\0';
        // So it's just "+" or "-"
        let strand = axt.q_strand.to_string();

        if let Some(psl) = Psl::from_align(
            &axt.q_name,
            q_size as u32,
            q_start,
            q_end,
            &axt.q_sym,
            &axt.t_name,
            t_size as u32,
            axt.t_start as i32,
            axt.t_end as i32,
            &axt.t_sym,
            &strand,
        ) {
            psl.write_to(&mut writer)?;
        }
    }

    Ok(())
}
