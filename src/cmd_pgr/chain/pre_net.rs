use anyhow::{bail, Result};
use clap::{Arg, ArgAction, ArgMatches, Command};
use pgr::libs::chaining::ChainReader;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};

pub fn make_subcommand() -> Command {
    Command::new("pre-net")
        .about("Remove chains that don't have a chance of being netted")
        .arg(Arg::new("input").required(true).help("Input chain file"))
        .arg(
            Arg::new("target_sizes")
                .required(true)
                .help("Target sizes file"),
        )
        .arg(
            Arg::new("query_sizes")
                .required(true)
                .help("Query sizes file"),
        )
        .arg(Arg::new("output").required(true).help("Output chain file"))
        .arg(
            Arg::new("dots")
                .long("dots")
                .value_parser(clap::value_parser!(usize))
                .help("Output a dot every so often"),
        )
        .arg(
            Arg::new("pad")
                .long("pad")
                .default_value("1")
                .value_parser(clap::value_parser!(u64))
                .help("Extra to pad around blocks to decrease trash"),
        )
        .arg(
            Arg::new("incl_hap")
                .long("incl-hap")
                .action(ArgAction::SetTrue)
                .help("Include query sequences name in the form *_hap*|*_alt*"),
        )
}

struct BitMap {
    size: u64,
    bits: Vec<u64>,
}

impl BitMap {
    fn new(size: u64) -> Self {
        let num_words = (size + 63) / 64;
        Self {
            size,
            bits: vec![0; num_words as usize],
        }
    }

    fn set_range(&mut self, start: u64, len: u64) {
        if len == 0 {
            return;
        }
        let end = (start + len).min(self.size);
        let start = start.min(self.size);
        if start >= end {
            return;
        }

        let start_word = (start / 64) as usize;
        let end_word = ((end - 1) / 64) as usize;

        let start_bit = start % 64;
        let end_bit = (end - 1) % 64;

        if start_word == end_word {
            let mask = (!0u64 << start_bit) & (!0u64 >> (63 - end_bit));
            self.bits[start_word] |= mask;
        } else {
            // First word
            self.bits[start_word] |= !0u64 << start_bit;

            // Middle words
            for i in (start_word + 1)..end_word {
                self.bits[i] = !0u64;
            }

            // Last word
            self.bits[end_word] |= !0u64 >> (63 - end_bit);
        }
    }

    fn is_fully_set(&self, start: u64, len: u64) -> bool {
        if len == 0 {
            return true;
        }
        let end = (start + len).min(self.size);
        let start = start.min(self.size);
        if start >= end {
            return true;
        }

        let start_word = (start / 64) as usize;
        let end_word = ((end - 1) / 64) as usize;

        let start_bit = start % 64;
        let end_bit = (end - 1) % 64;

        if start_word == end_word {
            let mask = (!0u64 << start_bit) & (!0u64 >> (63 - end_bit));
            return (self.bits[start_word] & mask) == mask;
        } else {
            // First word
            let mask1 = !0u64 << start_bit;
            if (self.bits[start_word] & mask1) != mask1 {
                return false;
            }

            // Middle words
            for i in (start_word + 1)..end_word {
                if self.bits[i] != !0u64 {
                    return false;
                }
            }

            // Last word
            let mask2 = !0u64 >> (63 - end_bit);
            if (self.bits[end_word] & mask2) != mask2 {
                return false;
            }
        }

        true
    }
}

fn load_sizes(path: &str) -> Result<HashMap<String, BitMap>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let name = parts[0].to_string();
            let size: u64 = parts[1].parse()?;
            map.insert(name, BitMap::new(size));
        }
    }
    Ok(map)
}

fn is_haplotype(name: &str) -> bool {
    name.ends_with("_hap")
        || name.ends_with("_alt")
        || name.contains("_hap")
        || name.contains("_alt")
    // UCSC implementation:
    // boolean haplotype(char *chrom)
    // {
    // return (stringIn("_hap", chrom) != NULL) || (stringIn("_alt", chrom) != NULL);
    // }
    // Wait, the C code used `stringIn` which matches substring, but `haplotype` usually implies suffix in some contexts.
    // However, `chainPreNet.c` calls `haplotype(chain->qName)`.
    // Let's stick to `contains("_hap") || contains("_alt")` to be safe and match `stringIn`.
}

pub fn execute(args: &ArgMatches) -> Result<()> {
    let input_path = args.get_one::<String>("input").unwrap();
    let target_sizes_path = args.get_one::<String>("target_sizes").unwrap();
    let query_sizes_path = args.get_one::<String>("query_sizes").unwrap();
    let output_path = args.get_one::<String>("output").unwrap();

    let dots = args.get_one::<usize>("dots").copied();
    let pad = args.get_one::<u64>("pad").copied().unwrap_or(1);
    let incl_hap = args.get_flag("incl_hap");

    let mut t_hash = load_sizes(target_sizes_path)?;
    let mut q_hash = load_sizes(query_sizes_path)?;

    let f = File::open(input_path)?;
    let mut reader = ChainReader::new(BufReader::new(f));

    let out_file = File::create(output_path)?;
    let mut writer = BufWriter::new(out_file);

    let mut last_score = f64::MAX;
    let mut count = 0;

    while let Some(res) = reader.next() {
        let chain = res?;

        // Check sort order
        if chain.header.score as f64 > last_score {
            bail!(
                "Input not sorted by score: {} > {}",
                chain.header.score,
                last_score
            );
        }
        last_score = chain.header.score as f64;

        if let Some(d) = dots {
            if count > 0 && count % d == 0 {
                eprint!(".");
            }
        }
        count += 1;

        // Filter haplotype
        if !incl_hap && is_haplotype(&chain.header.q_name) {
            continue;
        }

        let t_chrom = t_hash.get_mut(&chain.header.t_name).ok_or_else(|| {
            anyhow::anyhow!("Target sequence {} not found in sizes", chain.header.t_name)
        })?;

        // We need to access q_chrom as well. But Rust ownership prevents mutable borrowing both from the same HashMap if they were in the same map.
        // Luckily they are in different HashMaps (t_hash and q_hash).
        let q_chrom = q_hash.get_mut(&chain.header.q_name).ok_or_else(|| {
            anyhow::anyhow!("Query sequence {} not found in sizes", chain.header.q_name)
        })?;

        // Check used
        // Need to iterate blocks
        let blocks = chain.to_blocks();
        let mut any_open = false;

        for b in &blocks {
            // Check query
            if !q_chrom.is_fully_set(b.q_start, b.q_end - b.q_start) {
                any_open = true;
                break;
            }
            // Check target
            if !t_chrom.is_fully_set(b.t_start, b.t_end - b.t_start) {
                any_open = true;
                break;
            }
        }

        if any_open {
            chain.write(&mut writer)?;

            // Mark as used with pad
            for b in &blocks {
                // Apply pad
                // setWithPad(qChrom, b->qStart, b->qEnd);
                // setWithPad(tChrom, b->tStart, b->tEnd);

                // setWithPad logic:
                // s -= pad; if (s < 0) s = 0;
                // e += pad; if (e > size) e = size;
                // bitSetRange(bits, s, e-s);

                let q_s = b.q_start.saturating_sub(pad);
                let q_len = (b.q_end + pad).min(q_chrom.size) - q_s;
                q_chrom.set_range(q_s, q_len);

                let t_s = b.t_start.saturating_sub(pad);
                let t_len = (b.t_end + pad).min(t_chrom.size) - t_s;
                t_chrom.set_range(t_s, t_len);
            }
        }
    }

    if dots.is_some() {
        eprintln!();
    }

    Ok(())
}
