use clap::{Arg, ArgMatches, Command};
use pgr::libs::psl::Psl;
use std::io::{BufRead, Write};
use std::str::FromStr;

pub fn make_subcommand() -> Command {
    Command::new("tochain")
        .about("Convert PSL to Chain format")
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
                .help("Output Chain file")
                .default_value("stdout"),
        )
        .arg(
            Arg::new("fix_strand")
                .long("fix-strand")
                .short('f')
                .action(clap::ArgAction::SetTrue)
                .help("Fix '-' target strand by reverse complementing the record"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = args.get_one::<String>("input").unwrap();
    let output = args.get_one::<String>("output").unwrap();
    let fix_strand = args.get_flag("fix_strand");

    let reader = intspan::reader(input);
    let mut writer = intspan::writer(output);

    let mut chain_id = 1;

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let mut psl = match Psl::from_str(&line) {
            Ok(p) => p,
            Err(_) => {
                // Ignore lines that are not valid PSL (e.g. headers in some files)
                // Or maybe warn? UCSC pslFileOpen skips header lines if they look like header.
                // Here we assume standard PSL lines or skip errors.
                // Better to skip errors if it's just header.
                continue;
            }
        };

        // Handle negative target strand
        let t_strand_char = psl.strand.chars().nth(1).unwrap_or('+');
        if t_strand_char == '-' {
            if fix_strand {
                psl.rc();
            } else {
                // In strict mode we might abort, but for now maybe just warn or skip?
                // UCSC pslToChain aborts by default.
                // Let's abort to match behavior, or maybe just skip?
                // "errAbort" in C.
                anyhow::bail!("PSL record has '-' for target strand. Use --fix-strand to fix.");
            }
        }

        // Prepare Chain header fields
        let score = psl.score();
        let t_name = &psl.t_name;
        let t_size = psl.t_size;
        let t_start = psl.t_start;
        let t_end = psl.t_end;
        let q_name = &psl.q_name;
        let q_size = psl.q_size;
        let q_strand_char = psl.strand.chars().nth(0).unwrap_or('+');

        // Handle query strand for Chain format
        // Chain format: tStrand is always +, qStrand can be + or -
        // If qStrand is -, qStart/qEnd are relative to reverse end.
        let (q_start, q_end) = if q_strand_char == '-' {
            (q_size as i32 - psl.q_end, q_size as i32 - psl.q_start)
        } else {
            (psl.q_start, psl.q_end)
        };
        // Wait, if qStrand is '-', psl.qStart is coordinate on '+' strand.
        // Distance from end = size - coordinate.
        // If psl.qStart = 100, psl.qEnd = 200.
        // qStart (rev) = size - 200.
        // qEnd (rev) = size - 100.
        // So yes, (size - end, size - start).

        writeln!(
            writer,
            "chain {} {} {} + {} {} {} {} {} {} {} {}",
            score,
            t_name,
            t_size,
            t_start,
            t_end,
            q_name,
            q_size,
            q_strand_char,
            q_start,
            q_end,
            chain_id
        )?;

        // Write blocks
        for i in 0..psl.block_count as usize {
            let size = psl.block_sizes[i];
            write!(writer, "{}", size)?;

            if i < (psl.block_count as usize) - 1 {
                let dt = psl.t_starts[i + 1] - (psl.t_starts[i] + size);
                let dq = psl.q_starts[i + 1] - (psl.q_starts[i] + size);
                write!(writer, "\t{}\t{}", dt, dq)?;
            }
            writeln!(writer)?;
        }

        chain_id += 1;
    }

    Ok(())
}
