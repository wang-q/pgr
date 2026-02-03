use clap::*;
use pgr::libs::nt;
use pgr::libs::twobit::TwoBitFile;
use std::io::Write;

// Create clap subcommand arguments
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
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input 2bit file to process"),
        )
        .arg(
            Arg::new("ranges")
                .required(false)
                .index(2)
                .num_args(0..)
                .help("Ranges of interest"),
        )
        .arg(
            Arg::new("rgfile")
                .long("rgfile")
                .short('r')
                .num_args(1)
                .help("File of regions, one per line"),
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
    let infile = args.get_one::<String>("infile").unwrap();
    let output_path = args.get_one::<String>("outfile").unwrap();

    let mut ranges = if args.contains_id("ranges") {
        args.get_many::<String>("ranges")
            .unwrap()
            .cloned()
            .collect()
    } else {
        vec![]
    };

    if args.contains_id("rgfile") {
        let mut rgs = intspan::read_first_column(args.get_one::<String>("rgfile").unwrap());
        ranges.append(&mut rgs);
    }

    //----------------------------
    // Open files
    //----------------------------
    let mut tb = TwoBitFile::open(infile)?;
    let mut writer = intspan::writer(output_path);

    //----------------------------
    // Output
    //----------------------------
    for el in ranges.iter() {
        let rg = intspan::Range::from_str(el);
        let seq_id = rg.chr();

        // Check if sequence exists
        if !tb.sequence_offsets.contains_key(seq_id) {
            eprintln!("{} for [{}] not found in the 2bit file\n", seq_id, el);
            continue;
        }

        // Handle full sequence request (start=0 in intspan usually means just name)
        // intspan::Range::from_str("chr1") -> start=0, end=0
        let (start, end) = if *rg.start() == 0 {
            (None, None)
        } else {
            // Convert 1-based inclusive to 0-based half-open
            let s = (*rg.start() as usize).saturating_sub(1);
            let e = *rg.end() as usize;
            (Some(s), Some(e))
        };

        let mut seq = tb.read_sequence(seq_id, start, end, false)?;

        if rg.strand() == "-" {
            let rev_bytes: Vec<u8> = nt::rev_comp(seq.as_bytes()).collect();
            seq = String::from_utf8(rev_bytes).unwrap();
        }

        // Header construction
        let header = if *rg.start() == 0 {
            rg.to_string()
        } else {
            // Reconstruct range string with actual coordinates if needed,
            // but rg.to_string() usually gives what we want "chr:start-end"
            // If strand is -, intspan might output "chr(-):start-end"
            rg.to_string()
        };

        writeln!(writer, ">{}", header)?;
        writeln!(writer, "{}", seq)?;
    }

    Ok(())
}
