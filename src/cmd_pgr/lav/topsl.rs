use pgr::libs::lav::{LavReader, LavStanza, Block};
use pgr::libs::psl::Psl;
use clap::{Arg, ArgMatches, Command};
use std::fs::File;
use std::io::{self, BufReader, BufWriter};

pub fn make_subcommand() -> Command {
    Command::new("topsl")
        .about("Convert LAV to PSL format")
        .arg(
            Arg::new("input")
                .index(1)
                .help("Input LAV file (or stdin if not specified)")
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .help("Output PSL file (or stdout if not specified)")
                .num_args(1)
        )
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let input = matches.get_one::<String>("input");
    let output = matches.get_one::<String>("output");

    let reader: Box<dyn std::io::BufRead> = match input {
        Some(path) => Box::new(BufReader::new(File::open(path)?)),
        None => Box::new(BufReader::new(io::stdin())),
    };

    let mut writer: Box<dyn std::io::Write> = match output {
        Some(path) => Box::new(BufWriter::new(File::create(path)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };

    let mut lav_reader = LavReader::new(reader);
    
    // State
    let mut t_size = 0;
    let mut q_size = 0;
    let mut t_name = String::new();
    let mut q_name = String::new();
    let mut strand = String::from("+");

    while let Some(stanza) = lav_reader.next_stanza()? {
        match stanza {
            LavStanza::Sizes { t_size: t, q_size: q } => {
                t_size = t as u32;
                q_size = q as u32;
            }
            LavStanza::Header { t_name: t, q_name: q, is_rc } => {
                t_name = t;
                q_name = q;
                strand = if is_rc { "-".to_string() } else { "+".to_string() };
            }
            LavStanza::Alignment { blocks } => {
                if blocks.is_empty() { continue; }
                
                let psl = blocks_to_psl(&blocks, t_size, q_size, &t_name, &q_name, &strand);
                psl.write_to(&mut writer)?;
            }
            _ => {}
        }
    }

    Ok(())
}

fn blocks_to_psl(blocks: &[Block], t_size: u32, q_size: u32, t_name: &str, q_name: &str, strand: &str) -> Psl {
    let mut psl = Psl::new();
    psl.t_size = t_size;
    psl.q_size = q_size;
    psl.t_name = t_name.to_string();
    psl.q_name = q_name.to_string();
    psl.strand = strand.to_string();
    
    // Calculate overall range and stats
    let mut q_min = i64::MAX;
    let mut q_max = i64::MIN;
    let mut t_min = i64::MAX;
    let mut t_max = i64::MIN;

    for block in blocks {
        let len = (block.t_end - block.t_start) as u32;
        // UCSC lavToPsl calculation: match = (width * identity + 50)/100
        let match_cnt = (len * block.percent_id as u32 + 50) / 100;
        let mismatch_cnt = len - match_cnt;
        
        psl.match_count += match_cnt;
        psl.mismatch_count += mismatch_cnt;
        
        psl.block_count += 1;
        psl.block_sizes.push(len);
        psl.q_starts.push(block.q_start as u32);
        psl.t_starts.push(block.t_start as u32);

        if block.q_start < q_min { q_min = block.q_start; }
        if block.q_end > q_max { q_max = block.q_end; }
        if block.t_start < t_min { t_min = block.t_start; }
        if block.t_end > t_max { t_max = block.t_end; }
    }
    
    if !blocks.is_empty() {
        psl.q_start = q_min as i32;
        psl.q_end = q_max as i32;
        psl.t_start = t_min as i32;
        psl.t_end = t_max as i32;
    }

    // Gaps (inserts)
    for i in 0..blocks.len() - 1 {
        let curr = &blocks[i];
        let next = &blocks[i+1];
        
        // Assumption: blocks are sorted by T. LAV usually implies this.
        // If not, gap calculation might be weird (negative).
        // Let's assume non-negative gaps for now, or clamp to 0.
        
        let q_gap = next.q_start - curr.q_end;
        let t_gap = next.t_start - curr.t_end;
        
        if q_gap > 0 {
            psl.q_num_insert += 1;
            psl.q_base_insert += q_gap as i32;
        }
        
        if t_gap > 0 {
            psl.t_num_insert += 1;
            psl.t_base_insert += t_gap as i32;
        }
    }
    
    psl
}
