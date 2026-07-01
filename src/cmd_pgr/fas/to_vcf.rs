use anyhow::anyhow;
use clap::*;
use std::collections::{BTreeMap, HashSet};

use pgr::libs::alignment::{align_to_chr, get_subs, seq_intspan, vcf_alt_bases};
use pgr::libs::fmt::fas::next_fas_block;
use pgr::libs::fmt::vcf::{write_snp_row, write_vcf_header};

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
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("sizes")
                .long("sizes")
                .num_args(1)
                .help("Chrom sizes file with lines: <chr> <length>"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;
    let sizes_path = args
        .get_one::<String>("sizes")
        .map(|s| s.to_string())
        .unwrap_or_default();
    let sizes: BTreeMap<String, i32> = if !sizes_path.is_empty() {
        pgr::read_sizes(&sizes_path)?
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
                let mut reader = pgr::reader(infile)?;
                while let Ok(block) = next_fas_block(&mut reader) {
                    let chr = block
                        .entries
                        .first()
                        .ok_or_else(|| anyhow!("empty block entries in pre-scan"))?
                        .range()
                        .chr()
                        .to_string();
                    contigs.insert(chr);
                }
            }
        }
    }

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::reader(infile)?;

        while let Ok(block) = next_fas_block(&mut reader) {
            if !header_written {
                let contigs_ref = if sizes.is_empty() { None } else { Some(&sizes) };
                write_vcf_header(&mut writer, contigs_ref, &block.names)?;
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
            let target_entry = block
                .entries
                .get(target_entry_idx)
                .ok_or_else(|| anyhow!("missing target entry at index {}", target_entry_idx))?;
            let trange = target_entry.range().clone();
            let t_ints_seq = seq_intspan(target_entry.seq());

            let seq_count = seqs.len();
            let subs = get_subs(&seqs)?;

            for s in subs {
                let chr = trange.chr();
                let chr_pos = align_to_chr(&t_ints_seq, s.pos, trange.start, trange.strand())?;

                let pos_idx = (s.pos - 1) as usize;
                let ref_base = char::from(seqs[0][pos_idx]).to_ascii_uppercase();

                let alt_bases = vcf_alt_bases(&s);

                let sample_bases: Vec<u8> = seqs
                    .iter()
                    .take(seq_count)
                    .map(|seq| seq[pos_idx])
                    .collect();

                write_snp_row(
                    &mut writer,
                    chr,
                    chr_pos,
                    ref_base,
                    &alt_bases,
                    &sample_bases,
                )?;
            }
        }
    }

    Ok(())
}
