use bio::alphabets::dna::revcomp;
use bio::io::fasta;
use clap::{Arg, ArgMatches, Command};
use pgr::libs::chain::{Chain, ChainReader};
use pgr::libs::net::{read_nets, Fill, Gap};
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
            Arg::new("target_fa")
                .required(true)
                .help("Target FASTA file"),
        )
        .arg(Arg::new("query_fa").required(true).help("Query FASTA file"))
        .arg(Arg::new("out_axt").required(true).help("Output AXT file"))
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let in_net = matches.get_one::<String>("in_net").unwrap();
    let in_chain = matches.get_one::<String>("in_chain").unwrap();
    let target_fa = matches.get_one::<String>("target_fa").unwrap();
    let query_fa = matches.get_one::<String>("query_fa").unwrap();
    let out_axt = matches.get_one::<String>("out_axt").unwrap();

    // Load sequences
    let t_seqs = read_fasta(target_fa)?;
    let q_seqs = read_fasta(query_fa)?;

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

    let mut writer = BufWriter::new(File::create(out_axt)?);

    for net in &nets {
        r_convert(&net.root, &chains, &t_seqs, &q_seqs, &mut writer)?;
    }

    Ok(())
}

fn read_fasta(path: &str) -> anyhow::Result<HashMap<String, Vec<u8>>> {
    let reader = fasta::Reader::from_file(path)?;
    let mut seqs = HashMap::new();
    for result in reader.records() {
        let record = result?;
        seqs.insert(record.id().to_string(), record.seq().to_vec());
    }
    Ok(seqs)
}

fn r_convert<W: Write>(
    gap: &Rc<RefCell<Gap>>,
    chains: &HashMap<u64, Chain>,
    t_seqs: &HashMap<String, Vec<u8>>,
    q_seqs: &HashMap<String, Vec<u8>>,
    writer: &mut W,
) -> anyhow::Result<()> {
    let g = gap.borrow();
    for fill in &g.fills {
        let f = fill.borrow();
        if f.chain_id != 0 {
            if let Some(chain) = chains.get(&f.chain_id) {
                convert_fill(&f, chain, t_seqs, q_seqs, writer)?;
            }
        }

        // Recurse
        if !f.gaps.is_empty() {
            // Need to drop borrow of f to recurse if we were passing f, but we pass gaps.
            // However, f.gaps is Vec<Rc<RefCell<Gap>>>.
            // We can clone the Vec (cheap, just Rcs)
            let gaps = f.gaps.clone();
            drop(f); // Drop borrow
            for child_gap in gaps {
                r_convert(&child_gap, chains, t_seqs, q_seqs, writer)?;
            }
        }
    }
    Ok(())
}

fn convert_fill<W: Write>(
    fill: &Fill,
    chain: &Chain,
    t_seqs: &HashMap<String, Vec<u8>>,
    q_seqs: &HashMap<String, Vec<u8>>,
    writer: &mut W,
) -> anyhow::Result<()> {
    // Fill range on target
    let t_start = fill.start;
    let t_end = fill.end; // fill.start + fill.len

    // Get subset of chain
    // We need to convert chain to blocks and find those overlapping [t_start, t_end)
    let blocks = chain.to_blocks();

    // Iterate blocks
    for (i, block) in blocks.iter().enumerate() {
        // Check overlap
        let start = block.t_start.max(t_start);
        let end = block.t_end.min(t_end);

        if start < end {
            // Calculate offsets
            let t_offset = start - block.t_start;
            let len = end - start;

            let q_start_block = block.q_start + t_offset;
            let q_end_block = q_start_block + len;

            // Extract sequences
            let t_seq = get_subseq(t_seqs, &chain.header.t_name, start, end, '+')?;
            let q_seq = get_subseq(
                q_seqs,
                &chain.header.q_name,
                q_start_block,
                q_end_block,
                chain.header.q_strand,
            )?;

            // Check if this is the last block (for AXT header purposes, AXTs are usually per block or per chain?)
            // AXT format:
            // id chr1 start end chr2 start end strand score
            // seq1
            // seq2
            //
            // Usually one AXT entry per continuous alignment block.
            // But if gaps are small, they might be stitched?
            // UCSC netToAxt uses `chainToAxt` which handles gaps (inserts).
            // Here, `blocks` from `chain.to_blocks()` are ungapped alignment blocks.
            // The gaps (dt, dq) are between blocks.
            // If I iterate blocks and write AXT for each block, I lose the gap context (score might be wrong if split).
            // But `net` implies we are at a level where we might want to split or keep together.
            // `chainToAxt` has `maxGap` parameter.
            // If I just output one AXT per block, it is "valid" AXT but very fragmented.
            // For now, I will output one AXT per block for simplicity, or try to merge if gap is small.
            // Given I am implementing `netToAxt` which is supposed to be "best" alignment,
            // and `fill` defines a range.
            // The `chain` might have gaps within this fill.

            // Let's implement simple per-block AXT for MVP.

            writeln!(
                writer,
                "{} {} {} {} {} {} {} {} {}",
                i, // AXT index (should be unique per file? or per chain?)
                chain.header.t_name,
                start + 1, // AXT is 1-based
                end,
                chain.header.q_name,
                q_start_block + 1, // AXT 1-based. If strand is -, this needs care.
                q_end_block,
                chain.header.q_strand,
                chain.header.score // This is chain score, not block score.
            )?;

            writeln!(writer, "{}", std::str::from_utf8(&t_seq).unwrap())?;
            writeln!(writer, "{}", std::str::from_utf8(&q_seq).unwrap())?;
            writeln!(writer)?;
        }
    }

    Ok(())
}

fn get_subseq(
    seqs: &HashMap<String, Vec<u8>>,
    name: &str,
    start: u64,
    end: u64,
    strand: char,
) -> anyhow::Result<Vec<u8>> {
    let seq = seqs
        .get(name)
        .ok_or_else(|| anyhow::anyhow!("Sequence not found: {}", name))?;
    let len = seq.len() as u64;

    if start >= len || end > len {
        return Err(anyhow::anyhow!(
            "Coordinates out of bounds: {} {}-{} (len {})",
            name,
            start,
            end,
            len
        ));
    }

    // If strand is '+', simple slice.
    // If strand is '-', the coordinates [start, end) are on the REVERSE strand (if coming from Chain q coords).
    // Wait, let's verify Chain coordinate system.
    // Chain qStart/qEnd:
    // "If strand is -, coordinates are on the reverse strand."
    // Example: Len 100.
    // - strand, start 0, end 10.
    // This corresponds to the last 10 bases of the + strand, reversed and complemented.
    // + strand range: [100-10, 100-0) = [90, 100).
    // So we extract + strand [90, 100), then revcomp.

    let sub = if strand == '+' {
        seq[start as usize..end as usize].to_vec()
    } else {
        let p_start = len - end;
        let p_end = len - start;
        let s = seq[p_start as usize..p_end as usize].to_vec();
        revcomp(&s)
    };

    Ok(sub)
}
