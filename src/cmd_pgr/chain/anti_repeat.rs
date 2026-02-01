use clap::{Arg, ArgMatches, Command};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

use pgr::libs::chain::{read_chains, Chain, Block};

// Default scores from UCSC chainAntiRepeat.c
const DEFAULT_MIN_SCORE: usize = 5000;
const DEFAULT_NO_CHECK_SCORE: usize = 200000;

pub fn make_subcommand() -> Command {
    Command::new("anti-repeat")
        .about("Filter chains for repeats and degeneracy")
        .arg(
            Arg::new("target")
                .long("target")
                .short('t')
                .required(true)
                .help("Target sequence file (FASTA)"),
        )
        .arg(
            Arg::new("query")
                .long("query")
                .short('q')
                .required(true)
                .help("Query sequence file (FASTA)"),
        )
        .arg(
            Arg::new("input")
                .required(true)
                .help("Input chain file"),
        )
        .arg(
            Arg::new("output")
                .required(true)
                .help("Output chain file"),
        )
        .arg(
            Arg::new("min_score")
                .long("min-score")
                .default_value("5000")
                .value_parser(clap::value_parser!(usize))
                .help("Minimum score to pass"),
        )
        .arg(
            Arg::new("no_check_score")
                .long("no-check-score")
                .default_value("200000")
                .value_parser(clap::value_parser!(usize))
                .help("Score above which no checks are performed"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let target_path = args.get_one::<String>("target").unwrap();
    let query_path = args.get_one::<String>("query").unwrap();
    let input_path = args.get_one::<String>("input").unwrap();
    let output_path = args.get_one::<String>("output").unwrap();
    let min_score = *args.get_one::<usize>("min_score").unwrap();
    let no_check_score = *args.get_one::<usize>("no_check_score").unwrap();

    let target_seqs = load_fasta(target_path)?;
    let query_seqs = load_fasta(query_path)?;

    let mut reader = BufReader::new(File::open(input_path)?);
    let chains = read_chains(&mut reader)?; // Note: read_chains reads all chains into memory

    let mut writer = BufWriter::new(File::create(output_path)?);

    for chain in chains {
        if chain.header.score >= no_check_score as f64 {
            chain.write(&mut writer)?;
            continue;
        }

        if check_chain(&chain, &target_seqs, &query_seqs, min_score) {
            chain.write(&mut writer)?;
        }
    }

    Ok(())
}

fn load_fasta<P: AsRef<Path>>(path: P) -> anyhow::Result<HashMap<String, Vec<u8>>> {
    let reader = bio::io::fasta::Reader::from_file(path)?;
    let mut seqs = HashMap::new();
    for result in reader.records() {
        let record = result?;
        seqs.insert(record.id().to_string(), record.seq().to_vec());
    }
    Ok(seqs)
}

fn check_chain(
    chain: &Chain,
    t_seqs: &HashMap<String, Vec<u8>>,
    q_seqs: &HashMap<String, Vec<u8>>,
    min_score: usize,
) -> bool {
    let t_seq = match t_seqs.get(&chain.header.t_name) {
        Some(s) => s,
        None => return false,
    };
    let q_seq = match q_seqs.get(&chain.header.q_name) {
        Some(s) => s,
        None => return false,
    };

    let blocks = chain.to_blocks();
    
    // 1. Degeneracy Filter (Low complexity check)
    if !check_degeneracy(chain, &blocks, t_seq, q_seq, min_score) {
        return false;
    }

    // 2. Repeat Filter (Lowercase check)
    check_repeat(chain, &blocks, t_seq, q_seq, min_score)
}

fn get_slices<'a>(
    block: &Block, 
    t_seq: &'a [u8], 
    q_seq: &'a [u8], 
    q_strand: char, 
    q_size: u64
) -> Option<(&'a [u8], &'a [u8])> {
    if block.t_end as usize > t_seq.len() {
        return None;
    }
    let t_slice = &t_seq[block.t_start as usize..block.t_end as usize];

    let (q_start, q_end) = if q_strand == '+' {
        (block.q_start as usize, block.q_end as usize)
    } else {
        let q_len = q_size as usize;
        (q_len - block.q_end as usize, q_len - block.q_start as usize)
    };
    
    if q_end > q_seq.len() {
        return None;
    }
    let q_slice = &q_seq[q_start..q_end];
    
    Some((t_slice, q_slice))
}

fn check_degeneracy(
    chain: &Chain,
    blocks: &[Block],
    t_seq: &[u8],
    q_seq: &[u8],
    min_score: usize,
) -> bool {
    let mut counts = [0; 4]; // T, C, A, G
    let mut total_matches = 0;

    for block in blocks {
        if let Some((t_slice, q_slice)) = get_slices(block, t_seq, q_seq, chain.header.q_strand, chain.header.q_size) {
            for i in 0..t_slice.len() {
                let t_base = t_slice[i];
                let q_base_raw = if chain.header.q_strand == '+' {
                    q_slice[i]
                } else {
                    q_slice[q_slice.len() - 1 - i]
                };
                
                let t_val = nt_val(t_base);
                let mut q_val = nt_val(q_base_raw);

                if chain.header.q_strand == '-' && q_val >= 0 {
                    q_val = (q_val + 2) % 4;
                }
                
                if t_val >= 0 && t_val == q_val {
                    counts[t_val as usize] += 1;
                    total_matches += 1;
                }
            }
        }
    }

    if total_matches == 0 {
        return false;
    }

    // Sum of top 2
    let mut counts_vec = counts.to_vec();
    counts_vec.sort_unstable_by(|a, b| b.cmp(a)); // Descending
    let best2 = counts_vec[0] + counts_vec[1];

    let ok_best2 = 0.80;
    let observed_best2 = best2 as f64 / total_matches as f64;
    let over_ok = observed_best2 - ok_best2;
    let max_over_ok = 1.0 - ok_best2;

    if over_ok <= 0.0 {
        true
    } else {
        let adjust_factor = 1.01 - over_ok / max_over_ok;
        let adjusted_score = chain.header.score * adjust_factor;
        if adjusted_score < min_score as f64 {
            eprintln!("Chain {} filtered by degeneracy: score {} -> {}", chain.header.id, chain.header.score, adjusted_score);
            false
        } else {
            true
        }
    }
}

fn check_repeat(
    chain: &Chain,
    blocks: &[Block],
    t_seq: &[u8],
    q_seq: &[u8],
    min_score: usize,
) -> bool {
    let mut rep_count = 0;
    let mut total = 0;

    for block in blocks {
        if let Some((t_slice, q_slice)) = get_slices(block, t_seq, q_seq, chain.header.q_strand, chain.header.q_size) {
            for i in 0..t_slice.len() {
                let t_base = t_slice[i];
                let q_base = if chain.header.q_strand == '+' {
                    q_slice[i]
                } else {
                    q_slice[q_slice.len() - 1 - i]
                };
                
                if is_lower(t_base) || is_lower(q_base) {
                    rep_count += 1;
                }
            }
            total += t_slice.len();
        }
    }
    
    if total == 0 {
        return false;
    }

    let adjusted_score = chain.header.score * 2.0 * ((total - rep_count) as f64) / (total as f64);
    if adjusted_score < min_score as f64 {
        eprintln!("Chain {} filtered by repeat: score {} -> {} (rep {}/{})", chain.header.id, chain.header.score, adjusted_score, rep_count, total);
        false
    } else {
        true
    }
}

fn nt_val(base: u8) -> i8 {
    match base {
        b'T' | b't' => 0,
        b'C' | b'c' => 1,
        b'A' | b'a' => 2,
        b'G' | b'g' => 3,
        _ => -1,
    }
}

fn is_lower(base: u8) -> bool {
    base >= b'a' && base <= b'z'
}

