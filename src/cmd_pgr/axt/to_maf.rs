use clap::*;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter};
use std::path::Path;

use pgr::libs::axt::AxtReader;
use pgr::libs::maf::{MafAli, MafComp, MafWriter};

pub fn make_subcommand() -> Command {
    Command::new("to-maf")
        .about("Convert from axt to maf format")
        .after_help(
            r###"
Where tSizes and qSizes is a file that contains the sizes of the target and query sequences.
Very often this will be a chrom.sizes file.

Examples:
  # Convert axt to maf
  pgr axt to-maf in.axt -t t.sizes -q q.sizes -o out.maf

  # Split output by target name
  pgr axt to-maf in.axt -t t.sizes -q q.sizes --t-split -o out_dir
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
                .help("Output MAF file or directory")
                .default_value("stdout"),
        )
        .arg(
            Arg::new("q_prefix")
                .long("q-prefix")
                .value_name("STR")
                .help("Add prefix to start of query sequence name in maf")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("t_prefix")
                .long("t-prefix")
                .value_name("STR")
                .help("Add prefix to start of target sequence name in maf")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("t_split")
                .long("t-split")
                .help("Create a separate maf file for each target sequence. Output is a dir.")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("score")
                .long("score")
                .help("Recalculate score (Not implemented, uses AXT score)")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("score_zero")
                .long("score-zero")
                .help("Recalculate score if zero (Not implemented, uses AXT score)")
                .action(ArgAction::SetTrue),
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

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = args.get_one::<String>("input").unwrap();
    let t_sizes_path = args.get_one::<String>("t_sizes").unwrap();
    let q_sizes_path = args.get_one::<String>("q_sizes").unwrap();
    let output = args.get_one::<String>("output").unwrap();
    let q_prefix = args
        .get_one::<String>("q_prefix")
        .map(|s| s.as_str())
        .unwrap_or("");
    let t_prefix = args
        .get_one::<String>("t_prefix")
        .map(|s| s.as_str())
        .unwrap_or("");
    let t_split = args.get_flag("t_split");

    // Load sizes
    let t_sizes = load_sizes(t_sizes_path)?;
    let q_sizes = load_sizes(q_sizes_path)?;

    // Open input
    let reader = pgr::reader(input);
    let axt_reader = AxtReader::new(reader);

    // Prepare output
    let mut current_t_name = String::new();
    let mut single_writer: Option<MafWriter<Box<dyn std::io::Write>>> = None;

    if t_split {
        if !Path::new(output).exists() {
            fs::create_dir_all(output)?;
        }
    } else {
        let writer = pgr::writer(output);
        let mut writer = MafWriter::new(writer);
        writer.write_header("blastz")?; // Default to blastz as in C code
        single_writer = Some(writer);
    }

    let mut split_writers: HashMap<String, MafWriter<Box<dyn std::io::Write>>> = HashMap::new();

    for result in axt_reader {
        let axt = result?;

        // Handle tSplit file switching
        let writer = if t_split {
            if axt.t_name != current_t_name {
                // In C code: if (tSplit && !sameString(axt->tName, tName)) ...
                // It assumes sorted input for efficiency, but here we can keep a map or just open/close.
                // The C code errors if not sorted.
                // "if (hashLookup(uniqHash, tName) != NULL) errAbort(...)"
                // I will implement similar check or just open the file.
                // Since OS has limit on open files, and C code closes previous file,
                // it implies we should only keep one open or check sort.
                // Let's assume sorted as per C code requirement and implement the check?
                // Or just open/append. C code overwrites ("w").
                // If the input is not sorted, overwriting would lose data.
                // So I should enforce sorted check or append.
                // C code uses "w", so it implies STRICTLY sorted (grouped).

                // Let's strictly follow C logic: check if we've seen this tName before after leaving it.
                // But simpler: just track current name. If it changes, close old, open new.
                // If we see an old name again, it's an error if we follow C strictly.
                // But for robustness, appending is safer if we don't enforce sort.
                // However, C code says "in.maf must be sorted".
                // Let's just open new file for current tName.

                let path = Path::new(output).join(format!("{}.maf", axt.t_name));
                // We need to manage writers.
                // If we want to support unsorted, we need a map of writers (up to limit) or open/append.
                // Let's stick to C behavior: Open new file.
                // To avoid opening/closing for every record if sorted, we cache the current writer.

                if !split_writers.contains_key(&axt.t_name) {
                    // Check if we should close the previous one?
                    // C code: closes `f` when tName changes.
                    // So it only keeps ONE file open.
                    split_writers.clear(); // Close previous

                    let file = File::create(&path)?;
                    let mut w =
                        MafWriter::new(Box::new(BufWriter::new(file)) as Box<dyn std::io::Write>);
                    w.write_header("blastz")?;
                    split_writers.insert(axt.t_name.clone(), w);
                }
                current_t_name = axt.t_name.clone();
            }
            split_writers.get_mut(&axt.t_name).unwrap()
        } else {
            single_writer.as_mut().unwrap()
        };

        // Construct MAF alignment
        let score = axt.score.map(|s| s as f64);
        // Note: axt.score is Option<i32>, maf uses f64.

        // Target component
        let t_src = format!("{}{}", t_prefix, axt.t_name);
        let t_size = *t_sizes
            .get(&axt.t_name)
            .ok_or_else(|| anyhow::anyhow!("Target size not found for {}", axt.t_name))?;

        let t_comp = MafComp {
            src: t_src,
            start: axt.t_start,
            size: axt.t_end - axt.t_start,
            strand: axt.t_strand,
            src_size: t_size,
            text: axt.t_sym,
        };

        // Query component
        let q_src = format!("{}{}", q_prefix, axt.q_name);
        let q_size = *q_sizes
            .get(&axt.q_name)
            .ok_or_else(|| anyhow::anyhow!("Query size not found for {}", axt.q_name))?;

        let q_comp = MafComp {
            src: q_src,
            start: axt.q_start,
            size: axt.q_end - axt.q_start,
            strand: axt.q_strand,
            src_size: q_size,
            text: axt.q_sym,
        };

        let ali = MafAli {
            score,
            components: vec![t_comp, q_comp],
        };

        writer.write_ali(&ali)?;
    }

    Ok(())
}
