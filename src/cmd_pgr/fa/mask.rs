use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("mask")
        .about("Masks regions in FASTA file(s)")
        .after_help(
            r###"
This command masks regions in FASTA files based on a region file (BED/GFF/etc.).

Masking modes:
* Soft-masking (default): Convert to lowercase
* Hard-masking (--hard): Replace with N's

Input format (runlist.json):
{
    "seq1": "1-100,200-300",    # Mask positions 1-100 and 200-300
    "seq2": "50-150",           # Mask positions 50-150
    "seq3": "1-50,90-100,..."   # Multiple regions allowed
}

Notes:
* 1-based coordinates
* Inclusive ranges
* Sequences not in runlist remain unchanged
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* Invalid ranges are silently ignored

Examples:
1. Soft-mask regions:
   pgr fa mask input.fa regions.json -o output.fa

2. Hard-mask regions:
   pgr fa mask input.fa regions.json --hard -o output.fa

3. Process gzipped files:
   pgr fa mask input.fa.gz regions.json -o output.fa.gz

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input FASTA file to process",
        ))
        .arg(
            Arg::new("runlist")
                .required(true)
                .num_args(1)
                .index(2)
                .help("JSON file specifying regions to mask"),
        )
        .arg(
            Arg::new("hard")
                .long("hard")
                .action(ArgAction::SetTrue)
                .help("Hard-mask regions (replace with N's)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let mut fa_in = pgr::libs::fmt::fa::reader(args.get_one::<String>("infile").unwrap())?;

    let json = intspan::read_json(args.get_one::<String>("runlist").unwrap());
    let runlists = intspan::json2set(&json);

    let is_hard = args.get_flag("hard");

    let mut fa_out = pgr::libs::fmt::fa::writer(crate::cmd_pgr::args::get_outfile(args))?;

    //----------------------------
    // Process
    //----------------------------
    for result in fa_in.records() {
        let record = result?;
        let name = String::from_utf8(record.name().into())?;
        let seq = record.sequence();

        if !runlists.contains_key(&name) {
            fa_out.write_record(&record)?;
            continue;
        }

        // Get the regions to mask for this sequence
        let ints = runlists.get(&name).unwrap();
        let seq_str = String::from_utf8(seq[..].into())?;
        let seq_out = pgr::libs::fmt::fa::mask_sequence(&seq_str, ints, is_hard);

        //----------------------------
        // Output
        //----------------------------
        let record_out = pgr::libs::fmt::fa::new_record(&name, seq_out.as_bytes());
        fa_out.write_record(&record_out)?;
    }

    Ok(())
}
