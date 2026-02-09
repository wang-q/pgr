use clap::*;
use itertools::Itertools;
use std::collections::{BTreeMap, HashSet};
use std::io::Write;

use pgr::libs::alignment::{align_to_chr, get_subs, seq_intspan};
use pgr::libs::fas::next_fas_block;

pub fn make_subcommand() -> Command {
    Command::new("to-vcf")
        .about("Outputs VCF file (substitutions only)")
        .after_help(
            r###"
Outputs VCF file (substitutions only) from block FA files.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* Outputs substitutions only; ID/QUAL/FILTER/INFO are '.'
* CHROM/POS are derived from the target range; REF is the target base; ALT are non-REF bases
* Use `--sizes` to emit `##contig=<ID=...,length=...>` headers

Examples:
1. Output VCF from a block FASTA:
   pgr fas to-vcf tests/fasr/example.fas

2. Output VCF with contig headers:
   pgr fas to-vcf --sizes tests/fasr/S288c.chr.sizes tests/fasr/YDL184C.fas

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input block FA file(s) to process"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
        .arg(
            Arg::new("sizes")
                .long("sizes")
                .num_args(1)
                .help("Chrom sizes file with lines: <chr> <length>"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap());
    let sizes_path = args
        .get_one::<String>("sizes")
        .map(|s| s.to_string())
        .unwrap_or_default();
    let sizes: BTreeMap<String, i32> = if !sizes_path.is_empty() {
        intspan::read_sizes(&sizes_path)
    } else {
        BTreeMap::new()
    };

    let mut header_written = false;

    // Pre-scan contigs to emit ##contig lines when sizes provided and inputs are regular files
    let mut contigs: HashSet<String> = HashSet::new();
    if !sizes.is_empty() {
        if let Some(infiles) = args.get_many::<String>("infiles") {
            for infile in infiles {
                if infile.to_lowercase() == "stdin" {
                    continue;
                }
                let mut reader = pgr::reader(infile);
                while let Ok(block) = next_fas_block(&mut reader) {
                    let chr = block.entries.first().unwrap().range().chr().to_string();
                    contigs.insert(chr);
                }
            }
        }
    }

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::reader(infile);

        while let Ok(block) = next_fas_block(&mut reader) {
            if !header_written {
                if !sizes.is_empty() {
                    for (chr, len) in sizes.iter() {
                        let meta = format!("##contig=<ID={},length={}>\n", chr, len);
                        writer.write_all(meta.as_ref())?;
                    }
                }
                writer.write_all(b"##fileformat=VCFv4.2\n")?;
                writer.write_all(
                    b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
                )?;
                let mut header =
                    String::from("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
                for name in &block.names {
                    header.push('\t');
                    header.push_str(name);
                }
                header.push('\n');
                writer.write_all(header.as_ref())?;
                header_written = true;
            }

            let mut seqs: Vec<&[u8]> = vec![];
            for entry in &block.entries {
                seqs.push(entry.seq().as_ref());
            }
            if seqs.is_empty() {
                continue;
            }

            let target_entry_idx = 0usize;
            let trange = block.entries.get(target_entry_idx).unwrap().range().clone();
            let t_ints_seq = seq_intspan(block.entries.get(target_entry_idx).unwrap().seq());

            let seq_count = seqs.len();
            let subs = get_subs(&seqs)?;

            for s in subs {
                let chr = trange.chr();
                let chr_pos = align_to_chr(&t_ints_seq, s.pos, trange.start, trange.strand())?;

                let pos_idx = (s.pos - 1) as usize;
                let ref_base = char::from(seqs[0][pos_idx]).to_ascii_uppercase();

                let mut alt_bases: Vec<char> = vec![];
                for i in 0..seq_count {
                    let b = char::from(seqs[i][pos_idx]).to_ascii_uppercase();
                    if matches!(b, 'A' | 'C' | 'G' | 'T') && b != ref_base {
                        alt_bases.push(b);
                    }
                }
                alt_bases = alt_bases.into_iter().unique().collect();

                let alt_str = if alt_bases.is_empty() {
                    ".".to_string()
                } else {
                    alt_bases.iter().map(|c| c.to_string()).join(",")
                };

                let mut row = String::new();
                row.push_str(chr);
                row.push('\t');
                row.push_str(&chr_pos.to_string());
                row.push('\t');
                row.push_str(".\t");
                row.push(ref_base);
                row.push('\t');
                row.push_str(&alt_str);
                row.push_str("\t.\t.\t.\tGT");

                for i in 0..seq_count {
                    row.push('\t');
                    let b = char::from(seqs[i][pos_idx]).to_ascii_uppercase();
                    let gt = if !matches!(b, 'A' | 'C' | 'G' | 'T') {
                        ".".to_string()
                    } else if b == ref_base {
                        "0".to_string()
                    } else {
                        match alt_bases.iter().position(|&x| x == b) {
                            Some(idx) => (idx + 1).to_string(),
                            None => ".".to_string(),
                        }
                    };
                    row.push_str(&gt);
                }

                row.push('\n');
                writer.write_all(row.as_ref())?;
            }
        }
    }

    Ok(())
}
