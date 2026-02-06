use clap::{Arg, ArgMatches, Command};
use pgr::libs::chain::sub_matrix::SubMatrix;
use pgr::libs::chain::{Chain, ChainReader};
use pgr::libs::net::{read_nets, Fill, Gap};
use pgr::libs::twobit::TwoBitFile;
use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::rc::Rc;

pub fn make_subcommand() -> Command {
    Command::new("to-axt")
        .about("Convert net (and chain) to axt")
        .arg(Arg::new("in_net").required(true).help("Input net file"))
        .arg(Arg::new("in_chain").required(true).help("Input chain file"))
        .arg(
            Arg::new("target")
                .required(true)
                .help("Target 2bit file"),
        )
        .arg(Arg::new("query").required(true).help("Query 2bit file"))
        .arg(Arg::new("out_axt").required(true).help("Output AXT file"))
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let in_net = matches.get_one::<String>("in_net").unwrap();
    let in_chain = matches.get_one::<String>("in_chain").unwrap();
    let target = matches.get_one::<String>("target").unwrap();
    let query = matches.get_one::<String>("query").unwrap();
    let out_axt = matches.get_one::<String>("out_axt").unwrap();

    // Load sequences (TwoBitFile)
    let mut t_2bit = TwoBitFile::open(target)?;
    let mut q_2bit = TwoBitFile::open(query)?;

    // Load chains
    let mut chains = HashMap::new();
    let chain_reader = ChainReader::new(File::open(in_chain)?);
    for chain_res in chain_reader {
        let chain = chain_res?;
        chains.insert(chain.header.id, chain);
    }

    // Read nets
    let reader = BufReader::new(File::open(in_net)?);
    let nets = read_nets(reader)?;

    // Load scoring matrix (default HoxD55)
    let matrix = SubMatrix::hoxd55();

    let mut writer = BufWriter::new(File::create(out_axt)?);
    
    // Write headers from the first net (if any)
    if let Some(first_net) = nets.first() {
        for comment in &first_net.comments {
            writeln!(writer, "{}", comment)?;
        }
    }

    let mut counter = 0;

    for net in &nets {
        r_convert(
            &net.root,
            &chains,
            &mut t_2bit,
            &mut q_2bit,
            &matrix,
            &mut writer,
            &mut counter,
        )?;
    }

    Ok(())
}

fn r_convert<W: Write>(
    gap: &Rc<RefCell<Gap>>,
    chains: &HashMap<u64, Chain>,
    t_2bit: &mut TwoBitFile<BufReader<File>>,
    q_2bit: &mut TwoBitFile<BufReader<File>>,
    matrix: &SubMatrix,
    writer: &mut W,
    counter: &mut usize,
) -> anyhow::Result<()> {
    let g = gap.borrow();
    for fill in &g.fills {
        let f = fill.borrow();
        if f.chain_id != 0 {
            if let Some(chain) = chains.get(&f.chain_id) {
                convert_fill(&f, chain, chains, t_2bit, q_2bit, matrix, writer, counter)?;
            }
        } else {
            // If no chain, just recurse into gaps
            for gap_rc in &f.gaps {
                r_convert(gap_rc, chains, t_2bit, q_2bit, matrix, writer, counter)?;
            }
        }
    }
    Ok(())
}

fn convert_fill<W: Write>(
    fill: &Fill,
    chain: &Chain,
    chains: &HashMap<u64, Chain>,
    t_2bit: &mut TwoBitFile<BufReader<File>>,
    q_2bit: &mut TwoBitFile<BufReader<File>>,
    matrix: &SubMatrix,
    writer: &mut W,
    counter: &mut usize,
) -> anyhow::Result<()> {
    let mut cur = fill.start;

    // Iterate gaps to interleave segments and children
    for gap_rc in &fill.gaps {
        let (g_start, g_end, has_children, q_gap_size) = {
            let g = gap_rc.borrow();
            (g.start, g.end, !g.fills.is_empty(), g.o_end - g.o_start)
        };

        // Decision: Split or Merge?
        // Split if: has_children OR q_gap_size > 0 (double-sided gap)
        // Merge if: !has_children AND q_gap_size == 0 (single-sided gap / indel)
        
        let should_split = has_children || q_gap_size > 0;

        if should_split {
            // 1. Segment before gap
            if g_start > cur {
                convert_segment(
                    cur, g_start, chain, t_2bit, q_2bit, matrix, writer, counter,
                )?;
            }

            // 2. Recurse into gap
            r_convert(gap_rc, chains, t_2bit, q_2bit, matrix, writer, counter)?;

            // 3. Update cur to skip this gap
            cur = cur.max(g_end);
        } else {
            // Merge: We treat this gap as part of the alignment (indel).
            // We do NOT call convert_segment here, nor r_convert.
            // We effectively extend the current segment over this gap.
            // convert_segment will handle the gap by inserting dashes.
        }
    }

    // 3. Tail
    if cur < fill.end {
        convert_segment(
            cur, fill.end, chain, t_2bit, q_2bit, matrix, writer, counter,
        )?;
    }

    Ok(())
}

fn convert_segment<W: Write>(
    t_start: u64,
    t_end: u64,
    chain: &Chain,
    t_2bit: &mut TwoBitFile<BufReader<File>>,
    q_2bit: &mut TwoBitFile<BufReader<File>>,
    matrix: &SubMatrix,
    writer: &mut W,
    counter: &mut usize,
) -> anyhow::Result<()> {
    // Get subset of chain
    let blocks = chain.to_blocks();

    // Find blocks overlapping [t_start, t_end)
    let mut idx_start = None;
    let mut idx_end = None;

    for (i, block) in blocks.iter().enumerate() {
        if block.t_end > t_start && block.t_start < t_end {
            if idx_start.is_none() {
                idx_start = Some(i);
            }
            idx_end = Some(i);
        }
    }

    if idx_start.is_none() {
        return Ok(());
    }

    let idx_start = idx_start.unwrap();
    let idx_end = idx_end.unwrap();

    let mut t_seq_all = String::new();
    let mut q_seq_all = String::new();

    // Helper to read Q sequence considering strand
    let read_q = |start: u64, end: u64, q_2bit: &mut TwoBitFile<BufReader<File>>| -> anyhow::Result<String> {
        let (r_start, r_end) = if chain.header.q_strand == '-' {
            (
                chain.header.q_size - end,
                chain.header.q_size - start,
            )
        } else {
            (start, end)
        };
        let mut seq = q_2bit.read_sequence(
            &chain.header.q_name,
            Some(r_start as usize),
            Some(r_end as usize),
            false,
        )?;
        if chain.header.q_strand == '-' {
            let rev = pgr::libs::nt::rev_comp(seq.as_bytes()).collect();
            seq = String::from_utf8(rev).unwrap();
        }
        Ok(seq)
    };

    // Calculate initial q_start for the AXT record
    let q_start_out_base = if idx_start > 0 && t_start < blocks[idx_start].t_start {
        // We start in the gap before block[idx_start]
        let prev = &blocks[idx_start - 1];
        let gap_start_t = prev.t_end;
        // Check if we skipped dq
        if t_start <= gap_start_t {
            prev.q_end // We include dq
        } else {
            // We started inside dt. dq is on Q. dt is on T. They are independent.
            // If we start at t_start > gap_start_t, we are "late" in T.
            // Q is stuck at prev.q_end + dq (since dq happened).
            blocks[idx_start].q_start // which is prev.q_end + dq
        }
    } else {
        // We start in block[idx_start]
        let b = &blocks[idx_start];
        let offset = t_start.saturating_sub(b.t_start);
        b.q_start + offset
    };
    
    // Correct q_start logic:
    // If we start in a block, easy.
    // If we start in a gap:
    //   If we include dq, q_start is prev.q_end.
    //   If we exclude dq, q_start is prev.q_end + dq (== cur.q_start).
    
    // Let's refine inside the loop.

    for i in idx_start..=idx_end {
        let block = &blocks[i];
        
        // 1. Handle gap BEFORE this block (if i > idx_start, OR i == idx_start and we overlap the gap)
        if i > 0 {
            let prev = &blocks[i - 1];
            // Gap range on T: [prev.t_end, block.t_start)
            // Overlap with Fill: [max(gap_start, t_start), min(gap_end, t_end))
            let gap_start_t = prev.t_end;
            let gap_end_t = block.t_start;
            
            let overlap_start = gap_start_t.max(t_start);
            let overlap_end = gap_end_t.min(t_end);
            
            if overlap_start < overlap_end {
                // There is overlap with dt (T gap)
                // Append T bases
                let t_chunk = t_2bit.read_sequence(
                    &chain.header.t_name,
                    Some(overlap_start as usize),
                    Some(overlap_end as usize),
                    false
                )?;
                t_seq_all.push_str(&t_chunk);
                
                // Append Q dashes
                for _ in 0..(overlap_end - overlap_start) {
                    q_seq_all.push('-');
                }
            }
            
            // Handle dq (Q gap)
            // It occurs "between" blocks.
            // If Fill includes the boundary (prev.t_end), we include dq.
            // "Includes boundary": t_start <= prev.t_end < t_end
            if t_start <= gap_start_t && gap_start_t < t_end {
                let dq_len = block.q_start - prev.q_end;
                if dq_len > 0 {
                    let q_chunk = read_q(prev.q_end, block.q_start, q_2bit)?;
                    q_seq_all.push_str(&q_chunk);
                    
                    for _ in 0..dq_len {
                        t_seq_all.push('-');
                    }
                }
            }
        }
        
        // 2. Handle Block
        let start = block.t_start.max(t_start);
        let end = block.t_end.min(t_end);
        
        if start < end {
            let t_offset = start - block.t_start;
            let len = end - start;
            
            let t_chunk = t_2bit.read_sequence(
                &chain.header.t_name,
                Some(start as usize),
                Some(end as usize),
                false
            )?;
            t_seq_all.push_str(&t_chunk);
            
            let q_start_seg = block.q_start + t_offset;
            let q_end_seg = q_start_seg + len;
            let q_chunk = read_q(q_start_seg, q_end_seg, q_2bit)?;
            q_seq_all.push_str(&q_chunk);
        }
        
        // 3. Handle gap AFTER this block (only if this is the last processed block, check if Fill extends further)
        if i == idx_end {
             // Check if Fill extends beyond block.t_end
             if t_end > block.t_end {
                 // We might have a gap after this block that is partially covered
                 if i + 1 < blocks.len() {
                     let next = &blocks[i + 1];
                     let gap_start_t = block.t_end;
                     let gap_end_t = next.t_start;
                     
                     let overlap_start = gap_start_t.max(t_start);
                     let overlap_end = gap_end_t.min(t_end);
                     
                     if overlap_start < overlap_end {
                         let t_chunk = t_2bit.read_sequence(
                            &chain.header.t_name,
                            Some(overlap_start as usize),
                            Some(overlap_end as usize),
                            false
                        )?;
                        t_seq_all.push_str(&t_chunk);
                        for _ in 0..(overlap_end - overlap_start) {
                            q_seq_all.push('-');
                        }
                    }
                    
                    // Handle dq at block.t_end
                    // If Fill covers block.t_end
                    if t_start <= gap_start_t && gap_start_t < t_end {
                        let dq_len = next.q_start - block.q_end;
                        if dq_len > 0 {
                            let q_chunk = read_q(block.q_end, next.q_start, q_2bit)?;
                            q_seq_all.push_str(&q_chunk.to_ascii_uppercase());
                            for _ in 0..dq_len {
                                t_seq_all.push('-');
                            }
                        }
                    }
                 }
             }
        }
    }
    
    // Calculate final q_end based on q_seq content (bases only)
    let q_bases_count = q_seq_all.chars().filter(|c| *c != '-').count() as u64;
    let q_end_out = q_start_out_base + q_bases_count;

    // Calculate score
    let score = calculate_score(&t_seq_all, &q_seq_all, matrix);

    writeln!(
        writer,
        "{} {} {} {} {} {} {} {} {}",
        *counter,
        chain.header.t_name,
        t_start + 1, // AXT 1-based
        t_end,
        chain.header.q_name,
        q_start_out_base + 1, // AXT 1-based
        q_end_out,
        chain.header.q_strand,
        score
    )?;
    *counter += 1;

    writeln!(writer, "{}", t_seq_all)?;
    writeln!(writer, "{}", q_seq_all)?;
    writeln!(writer)?;

    Ok(())
}

fn calculate_score(t_seq: &str, q_seq: &str, matrix: &SubMatrix) -> i32 {
    let mut score = 0;
    let t_chars: Vec<char> = t_seq.chars().collect();
    let q_chars: Vec<char> = q_seq.chars().collect();
    let len = t_chars.len();
    
    let mut in_gap_t = false;
    let mut in_gap_q = false;
    
    for i in 0..len {
        let t = t_chars[i];
        let q = q_chars[i];
        
        if t == '-' {
            // Gap in T (insertion in Q)
            if !in_gap_t {
                score -= matrix.gap_open;
                in_gap_t = true;
            }
            score -= matrix.gap_extend;
            in_gap_q = false; 
        } else if q == '-' {
            // Gap in Q (deletion in Q)
            if !in_gap_q {
                score -= matrix.gap_open;
                in_gap_q = true;
            }
            score -= matrix.gap_extend;
            in_gap_t = false;
        } else {
            // Match/Mismatch
            score += matrix.get_score(t, q);
            in_gap_t = false;
            in_gap_q = false;
        }
    }
    score
}
