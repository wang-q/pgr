use pgr::libs::chaining::{calc_block_score, chain_blocks, ChainableBlock, ScoreContext};
use pgr::libs::chaining::GapCalc;
use pgr::libs::chaining::SubMatrix;
use pgr::libs::psl::Psl;
use pgr::libs::twobit::TwoBitFile;
use anyhow::Result;
use clap::{Arg, ArgMatches, Command};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::io::{BufRead, Write};
use std::str::FromStr;

pub fn make_subcommand() -> Command {
    Command::new("psl")
        .about("Chain PSL alignments")
        .after_help(
            r###"
Processing:
  1. Group PSL blocks by target/query sequence and strand.
  2. Build a KD-tree (k-dimensional tree) for efficient predecessor search.
     - In this context, it's a 2D tree indexing blocks by (query_start, target_start).
     - It allows fast range queries to find candidate predecessor blocks that are "before" the current block in both query and target coordinates.
  3. Connect blocks into chains using dynamic programming:
     - Score = BlockScore + Max(PredecessorScore - GapCost).
     - GapCost depends on distance (linear: loose/medium).
     - Overlaps are trimmed by finding the optimal cut point based on exact sequence scores.
  4. Filter chains by minimum score.

Examples:
  # Chain PSL file
  pgr chaining psl t.2bit q.2bit in.psl -o out.chain
"###,
        )
        .arg(
            Arg::new("target")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Path to the target genome 2bit file"),
        )
        .arg(
            Arg::new("query")
                .required(true)
                .num_args(1)
                .index(2)
                .help("Path to the query genome 2bit file"),
        )
        .arg(
            Arg::new("psl")
                .required(true)
                .num_args(1)
                .index(3)
                .help("Path to the PSL file"),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("FILE")
                .help("Output Chain file")
                .default_value("-"),
        )
        .arg(
            Arg::new("linear_gap")
                .long("linear-gap")
                .default_value("medium")
                .value_parser(["loose", "medium"])
                .help("Linear gap cost type"),
        )
        .arg(
            Arg::new("min_score")
                .long("min-score")
                .default_value("0")
                .value_parser(clap::value_parser!(f64))
                .help("Minimum score of chain"),
        )
        .arg(
            Arg::new("gap_open")
                .long("gap-open")
                .value_parser(clap::value_parser!(i32))
                .help("Gap open cost (overrides --linear-gap)"),
        )
        .arg(
            Arg::new("gap_extend")
                .long("gap-extend")
                .value_parser(clap::value_parser!(i32))
                .help("Gap extension cost (overrides --linear-gap)"),
        )
        .arg(
            Arg::new("score_scheme")
                .long("score-scheme")
                .value_name("FILE")
                .help("Score scheme file (LASTZ format)"),
        )
}

pub fn execute(args: &ArgMatches) -> Result<()> {
    let input = args.get_one::<String>("psl").unwrap();
    let output = args.get_one::<String>("output").unwrap();
    let linear_gap = args.get_one::<String>("linear_gap").unwrap();
    let min_score = *args.get_one::<f64>("min_score").unwrap();
    let target_2bit_path = args.get_one::<String>("target");
    let query_2bit_path = args.get_one::<String>("query");
    let score_scheme_path = args.get_one::<String>("score_scheme");

    let reader = intspan::reader(input);
    let mut writer = intspan::writer(output);

    let mut t_2bit = if let Some(path) = target_2bit_path {
        Some(TwoBitFile::open(path)?)
    } else {
        None
    };

    let mut q_2bit = if let Some(path) = query_2bit_path {
        Some(TwoBitFile::open(path)?)
    } else {
        None
    };

    let score_matrix = if let Some(path) = score_scheme_path {
                SubMatrix::from_file(path)?
            } else {
                SubMatrix::default()
    };
    
    let mut score_context = if t_2bit.is_some() && q_2bit.is_some() {
        Some(ScoreContext {
            t_2bit: t_2bit.as_mut().unwrap(),
            q_2bit: q_2bit.as_mut().unwrap(),
            matrix: &score_matrix,
        })
    } else {
        None
    };

    let gap_open = args.get_one::<i32>("gap_open");
    let gap_extend = args.get_one::<i32>("gap_extend");

    let gap_calc = if let (Some(&open), Some(&extend)) = (gap_open, gap_extend) {
        GapCalc::affine(open, extend)
    } else {
        match linear_gap.as_str() {
            "loose" => GapCalc::loose(),
            "medium" => GapCalc::medium(),
            _ => GapCalc::medium(),
        }
    };

    // Group blocks by (t_name, q_name, q_strand)
    let mut groups: HashMap<(String, String, char), GroupData> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let psl = match Psl::from_str(&line) {
            Ok(p) => p,
            Err(_) => continue, // Skip invalid lines (e.g. headers)
        };

        let t_name = psl.t_name.clone();
        let q_name = psl.q_name.clone();
        let q_strand = psl.strand.chars().nth(0).unwrap_or('+');

        let key = (t_name.clone(), q_name.clone(), q_strand);
        let entry = groups.entry(key).or_insert_with(|| GroupData {
            t_size: psl.t_size,
            q_size: psl.q_size,
            blocks: Vec::new(),
        });

        if psl.strand.len() > 1 && psl.strand.chars().nth(1) == Some('-') {
            eprintln!(
                "Warning: Skipping PSL record with negative target strand: {} {} {}",
                psl.q_name, psl.strand, psl.t_name
            );
            continue;
        }

        for i in 0..psl.block_count as usize {
            let size = psl.block_sizes[i] as u64;
            let t_start = psl.t_starts[i] as u64;
            let t_end = t_start + size;

            let (q_start, q_end) = {
                let s = psl.q_starts[i] as u64;
                (s, s + size)
            };
            
            let mut block = ChainableBlock {
                t_start,
                t_end,
                q_start,
                q_end,
                score: size as f64 * 100.0,
            };

            if let Some(ctx) = &mut score_context {
                if let Some(exact) = calc_block_score(
                    &block,
                    ctx,
                    &q_name,
                    &t_name,
                    psl.q_size as u64,
                    q_strand,
                ) {
                    block.score = exact;
                }
            }

            entry.blocks.push(block);
        }
    }

    // Process groups
    let mut all_chains = Vec::new();
    let mut chain_id_counter = 1;

    for ((t_name, q_name, q_strand), mut data) in groups {
        if data.blocks.is_empty() {
            continue;
        }

        data.blocks.sort_by(|a, b| a.t_start.cmp(&b.t_start));

        if std::env::var("PGR_DEBUG").is_ok() {
            eprintln!("Group: {} {} {}", t_name, q_name, q_strand);
            for b in &data.blocks {
                eprintln!("Block: T {}-{} Q {}-{} Score {}", b.t_start, b.t_end, b.q_start, b.q_end, b.score);
            }
        }

        let chains = chain_blocks(
            &data.blocks,
            &gap_calc,
            &mut score_context,
            &q_name,
            data.q_size as u64,
            q_strand,
            &t_name,
            data.t_size as u64,
            &mut chain_id_counter,
        );
        all_chains.extend(chains);
    }

    all_chains.sort_by(|a, b| b.header.score.partial_cmp(&a.header.score).unwrap_or(Ordering::Equal));

    for chain in all_chains {
        if chain.header.score < min_score {
            continue;
        }

        write!(
            writer,
            "chain {:.0} {} {} {} {} {} {} {} {} {} {} {}\n",
            chain.header.score,
            chain.header.t_name,
            chain.header.t_size,
            chain.header.t_strand,
            chain.header.t_start,
            chain.header.t_end,
            chain.header.q_name,
            chain.header.q_size,
            chain.header.q_strand,
            chain.header.q_start,
            chain.header.q_end,
            chain.header.id
        )?;

        for (i, d) in chain.data.iter().enumerate() {
            if i == chain.data.len() - 1 {
                write!(writer, "{}\n", d.size)?;
            } else {
                write!(writer, "{}\t{}\t{}\n", d.size, d.dt, d.dq)?;
            }
        }
        write!(writer, "\n")?;
    }

    Ok(())
}

struct GroupData {
    t_size: u32,
    q_size: u32,
    blocks: Vec<ChainableBlock>,
}
