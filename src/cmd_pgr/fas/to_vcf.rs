use anyhow::{anyhow, Context};
use clap::{Arg, ArgMatches, Command};
use std::collections::BTreeMap;
use std::io::Write;

use pgr::libs::alignment::{align_to_chr, get_subs, seq_intspan, vcf_alt_bases};
use pgr::libs::fmt::fas::iter_fas_blocks;
use pgr::libs::fmt::vcf::{write_snp_row, write_vcf_header};
/// Build the clap subcommand for to-vcf.
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
/// Execute the to-vcf command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;
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

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader =
            pgr::reader(infile).with_context(|| format!("open reader for {}", infile))?;

        let mut block_idx = 0usize;
        for block_result in iter_fas_blocks(&mut reader) {
            let block = block_result
                .with_context(|| format!("read block {} from {}", block_idx, infile))?;
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
                block_idx += 1;
                continue;
            }

            let target_entry_idx = 0usize;
            let target_entry = block.entries.get(target_entry_idx).ok_or_else(|| {
                anyhow!(
                    "missing target entry at index {} in block {}",
                    target_entry_idx,
                    block_idx
                )
            })?;
            let trange = target_entry.range();
            let t_ints_seq = seq_intspan(target_entry.seq());

            let seq_count = seqs.len();
            let subs = get_subs(&seqs)?;

            for s in subs {
                let chr = trange.chr();
                let chr_pos = align_to_chr(&t_ints_seq, s.pos, trange.start, trange.strand())
                    .with_context(|| {
                        format!("align_to_chr at pos {} in block {}", s.pos, block_idx)
                    })?;

                let pos_idx = usize::try_from(s.pos).map_err(|_| {
                    anyhow!("invalid substitution pos {} in block {}", s.pos, block_idx)
                })?;
                let pos_idx = pos_idx.checked_sub(1).ok_or_else(|| {
                    anyhow!("invalid substitution pos {} in block {}", s.pos, block_idx)
                })?;
                if pos_idx >= seqs[0].len() {
                    anyhow::bail!(
                        "substitution pos {} out of range (seq len {}) in block {}",
                        s.pos,
                        seqs[0].len(),
                        block_idx
                    );
                }
                let ref_base = char::from(seqs[0][pos_idx]).to_ascii_uppercase();

                let alt_bases = vcf_alt_bases(&s);

                let sample_bases: Vec<u8> = seqs
                    .iter()
                    .take(seq_count)
                    .map(|seq| seq.get(pos_idx).copied().unwrap_or(b'-'))
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
            block_idx += 1;
        }
    }

    writer.flush()?;
    Ok(())
}
