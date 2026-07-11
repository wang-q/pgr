use anyhow::Context;
use clap::{ArgMatches, Command};
use pgr::libs::pbit::decompressor::Decompressor;
use std::io::Write;

/// Build the clap subcommand for range.
pub fn make_subcommand() -> Command {
    Command::new("range")
        .about("Extracts sequence regions by coordinates from all samples")
        .after_help(
            r###"
This command extracts sequence regions from a pbit archive using genomic
coordinates. For each region, it iterates over ALL samples that contain
the contig and writes one FASTA entry per sample (getctg semantics).

Range format:
    seq_name(strand):start-end

* seq_name: Required, contig identifier
* strand: Optional, + (default) or -
* start-end: Required, 1-based coordinates

Notes:
* All coordinates (<start> and <end>) are based on the positive strand,
  regardless of the specified strand
* Coordinates are 1-based inclusive
* pbit files are binary and require random access (seeking)
* Does not support stdin or gzipped inputs

Examples:
1. Extract ranges from command line:
   pgr pbit range input.pbit "chr1:1-1000" "chr1(+):90-150"

2. Extract ranges from file:
   pgr pbit range input.pbit -r ranges.txt

3. Extract with negative strand:
   pgr pbit range input.pbit "chr1(-):1-1000"
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input pbit file to process",
        ))
        .arg(crate::cmd_pgr::args::ranges_arg())
        .arg(crate::cmd_pgr::args::rgfile_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the range command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args
        .get_one::<String>("infile")
        .context("missing required argument: infile")?;
    let output_path = crate::cmd_pgr::args::get_outfile(args);

    let ranges = crate::cmd_pgr::args::collect_ranges(args)?;

    let mut dec = Decompressor::open(infile)
        .with_context(|| format!("Failed to open pbit file {}", infile))?;
    let mut writer = pgr::libs::io::writer(output_path)
        .with_context(|| format!("Failed to open writer for {}", output_path))?;

    for el in ranges.iter() {
        let rg = intspan::Range::from_str(el);
        let contig = rg.chr();

        // A range without ':' is a full-contig request (e.g. "chr1");
        // intspan returns start=0/end=0 and is_valid=false for these, so
        // bypass validation. Anything with ':' must parse as a valid range.
        let is_full_contig = !el.contains(':');
        if !is_full_contig && !rg.is_valid() {
            log::warn!("invalid range format: {}", el);
            continue;
        }

        // Check if contig exists in any sample.
        if !dec.contains_contig(contig) {
            log::warn!("{} for [{}] not found in any sample", contig, el);
            continue;
        }

        // Handle full contig request (start=0 and end=0 means just the name).
        let (start, end) = if is_full_contig {
            (None, None)
        } else {
            let start_val = *rg.start();
            let end_val = *rg.end();
            anyhow::ensure!(
                start_val > 0 && end_val > 0,
                "range coordinates must be positive: {}",
                el
            );
            // Convert 1-based inclusive to 0-based half-open.
            let s = (start_val as usize).saturating_sub(1);
            let e = end_val as usize;
            (Some(s), Some(e))
        };

        let strand = if rg.strand() == "-" { "-" } else { "+" };
        dec.get_contig(contig, start, end, strand, &mut writer)?;
    }

    writer.flush()?;
    Ok(())
}
