use clap::*;
use std::io::Write;
use pgr::libs::twobit::TwoBitFile;
use std::ops::Range;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("masked")
        .about("Identify masked regions in 2bit file(s)")
        .after_help(
            r###"
This command identifies masked regions in one or more 2bit files. Masked regions can be:
- Soft-masked regions (lowercase in FASTA, stored as mask blocks in 2bit)
- Hard-masked regions (N/n in FASTA, stored as N blocks in 2bit)

The output is a list of regions in the format:
    seq_name:start-end        # For regions spanning multiple positions
    seq_name:position        # For single positions

Notes:
* Coordinates are 1-based, inclusive
* Adjacent masked positions are merged into a single region

Examples:
1. Identify masked regions (soft and hard):
   pgr 2bit masked input.2bit -o masked_regions.txt

2. Identify only N/n gap regions (hard-masked):
   pgr 2bit masked input.2bit --gap -o gap_regions.txt

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input 2bit file(s) to process"),
        )
        .arg(
            Arg::new("gap")
                .long("gap")
                .short('g')
                .action(ArgAction::SetTrue)
                .help("Only identify regions of N/n (gaps)"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let is_gap = args.get_flag("gap");
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    //----------------------------
    // Ops
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut tb = TwoBitFile::open(infile)?;
        let names = tb.get_sequence_names();

        for name in names {
            let (n_blocks, mask_blocks) = tb.get_sequence_blocks(&name)?;
            
            let mut blocks = Vec::new();
            blocks.extend(n_blocks.0.into_iter());
            if !is_gap {
                blocks.extend(mask_blocks.0.into_iter());
            }

            if blocks.is_empty() {
                continue;
            }

            // Sort by start position
            blocks.sort_by(|a, b| a.start.cmp(&b.start));

            // Merge adjacent or overlapping blocks
            let mut merged: Vec<Range<usize>> = Vec::new();
            if let Some(first) = blocks.first() {
                merged.push(first.clone());
            }

            for block in blocks.iter().skip(1) {
                let last = merged.last_mut().unwrap();
                // Check for overlap or adjacency
                // block.start <= last.end handles both:
                // Overlap: [0, 5) and [3, 8) -> 3 < 5
                // Adjacency: [0, 5) and [5, 10) -> 5 == 5
                if block.start <= last.end {
                    last.end = last.end.max(block.end);
                } else {
                    merged.push(block.clone());
                }
            }

            // Write output
            for block in merged {
                // Convert 0-based half-open [start, end) to 1-based inclusive
                let start_1based = block.start + 1;
                let end_1based = block.end;
                
                if start_1based == end_1based {
                    writer.write_all(format!("{}:{}\n", name, start_1based).as_bytes())?;
                } else {
                    writer.write_all(format!("{}:{}-{}\n", name, start_1based, end_1based).as_bytes())?;
                }
            }
        }
    }

    Ok(())
}
