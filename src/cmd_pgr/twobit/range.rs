use anyhow::Context;
use clap::{ArgMatches, Command};
use pgr::libs::fmt::twobit::TwoBitFile;
use pgr::libs::nt;
use std::io::Write;

/// Build the clap subcommand for range.
pub fn make_subcommand() -> Command {
    Command::new("range")
        .about("Extracts sequence regions by coordinates")
        .after_help(
            r###"
This command extracts sequence regions from 2bit files using genomic coordinates.

Range format:
    seq_name(strand):start-end

* seq_name: Required, sequence identifier
* strand: Optional, + (default) or -
* start-end: Required, 1-based coordinates

Notes:
* All coordinates (<start> and <end>) are based on the positive strand, regardless of the specified strand
* 2bit files support efficient random access, so no cache is needed
* 2bit files are binary and require random access (seeking)
* Does not support stdin or gzipped inputs

Examples:
1. Extract ranges from command line:
   pgr 2bit range input.2bit "chr1:1-1000" "chr1(+):90-150"

2. Extract ranges from file:
   pgr 2bit range input.2bit -r ranges.txt

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input 2bit file to process",
        ))
        .arg(crate::cmd_pgr::args::ranges_arg())
        .arg(crate::cmd_pgr::args::rgfile_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the range command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args.get_one::<String>("infile").unwrap();
    let output_path = crate::cmd_pgr::args::get_outfile(args);

    let ranges = crate::cmd_pgr::args::collect_ranges(args)?;

    // Open files
    let mut tb =
        TwoBitFile::open(infile).with_context(|| format!("Failed to open 2bit file {}", infile))?;
    let mut writer = pgr::writer(output_path)
        .with_context(|| format!("Failed to open writer for {}", output_path))?;

    for el in ranges.iter() {
        let rg = intspan::Range::from_str(el);
        let seq_id = rg.chr();

        // Check if sequence exists
        if !tb.sequence_offsets.contains_key(seq_id) {
            log::warn!("{} for [{}] not found in the 2bit file", seq_id, el);
            continue;
        }

        // Handle full sequence request (start=0 in intspan usually means just name)
        // intspan::Range::from_str("chr1") -> start=0, end=0
        let (start, end) = if *rg.start() == 0 {
            (None, None)
        } else {
            anyhow::ensure!(rg.is_valid(), "invalid range: {}", el);
            // Convert 1-based inclusive to 0-based half-open
            let s = (*rg.start() as usize).saturating_sub(1);
            let e = *rg.end() as usize;
            (Some(s), Some(e))
        };

        let mut seq = tb.read_sequence(seq_id, start, end, false)?;

        if rg.strand() == "-" {
            let rev_bytes: Vec<u8> = nt::rev_comp(seq.as_bytes()).collect();
            seq = String::from_utf8(rev_bytes)
                .map_err(|e| anyhow::anyhow!("invalid utf8 in rev_comp: {}", e))?;
        }

        // Header construction
        let header = rg.to_string();

        writeln!(writer, ">{}", header)?;
        writeln!(writer, "{}", seq)?;
    }

    writer.flush()?;
    Ok(())
}
