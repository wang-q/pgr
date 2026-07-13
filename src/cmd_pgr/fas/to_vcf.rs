use anyhow::Context;
use clap::{Arg, ArgMatches, Command};
use std::collections::BTreeMap;
use std::io::Write;

use pgr::libs::fmt::fas::iter_fas_blocks;
use pgr::libs::fmt::vcf::write_vcf_header;

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
   pgr fas to-vcf tests/fas/example.fas

2. Output VCF with contig headers:
   pgr fas to-vcf --sizes tests/fas_vcf/S288c.chr.sizes tests/fas_vcf/YDL184C.fas

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
    let mut header_names: Option<Vec<String>> = None;

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader =
            pgr::reader(infile).with_context(|| format!("open reader for {}", infile))?;

        for (block_idx, block_result) in iter_fas_blocks(&mut reader).enumerate() {
            let block = block_result
                .with_context(|| format!("read block {} from {}", block_idx, infile))?;
            if !header_written {
                let contigs_ref = if sizes.is_empty() { None } else { Some(&sizes) };
                write_vcf_header(&mut writer, contigs_ref, &block.names)?;
                header_names = Some(block.names.clone());
                header_written = true;
            } else if let Some(ref expected) = header_names {
                if block.names != *expected {
                    anyhow::bail!(
                        "block {} from {} has inconsistent samples: expected {:?}, got {:?}; VCF requires the same species in the same order across all blocks",
                        block_idx,
                        infile,
                        expected,
                        block.names
                    );
                }
            }

            pgr::libs::fmt::fas::write_vcf_block(&block, block_idx, &mut writer)?;
        }
    }

    writer.flush()?;
    Ok(())
}
