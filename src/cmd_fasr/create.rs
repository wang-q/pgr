use clap::*;
use std::io::BufRead;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("create")
        .about("Create block FA files from links of ranges")
        .after_help(
            r###"
This subcommand creates block FA files from links of ranges, typically generated by `linkr`.

Input files can be gzipped. If the input file is 'stdin', data is read from standard input.

Note:
- The reference genome(s) must be provided as a multi-sequence FA file, can be bgzipped.
- Two styles of FA headers are supported:
  - `>chr` for single-genome self-alignments.
  - `>name.chr` for multiple genomes.
- Requires `samtools` to be installed and available in $PATH.

Examples:
1. Create block FA files for a single genome:
   fasr create tests/fasr/genome.fa tests/fasr/I.connect.tsv

2. Create block FA files for a specific species:
   fasr create tests/fasr/genome.fa tests/fasr/I.connect.tsv --name S288c

3. Create block FA files for multiple genomes:
   fasr create tests/fasr/genomes.fa tests/fasr/I.connect.tsv --multi

4. Output results to a file:
   fasr create tests/fasr/genome.fa tests/fasr/I.connect.tsv -o output.fas

"###,
        )
        .arg(
            Arg::new("genome.fa")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Path to the reference genome FA file"),
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(2)
                .help("Input file(s) containing links of ranges"),
        )
        .arg(
            Arg::new("name")
                .long("name")
                .num_args(1)
                .help("Set a species name for ranges. No effects if --multi"),
        )
        .arg(
            Arg::new("multi")
                .long("multi")
                .action(ArgAction::SetTrue)
                .help("Indicates the reference genome contains multiple genomes"),
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
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let opt_genome = args.get_one::<String>("genome.fa").unwrap();
    let opt_name = &args
        .get_one::<String>("name")
        .map(|s| s.as_str())
        .unwrap_or("")
        .to_string();
    let is_multi = args.get_flag("multi");

    //----------------------------
    // Ops
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = intspan::reader(infile);
        for line in reader.lines().map_while(Result::ok) {
            let parts: Vec<&str> = line.split('\t').collect();

            for part in &parts {
                let mut range = intspan::Range::from_str(part);
                if !range.is_valid() {
                    continue;
                }

                // Set the species name if provided
                if !opt_name.is_empty() {
                    *range.name_mut() = opt_name.to_string();
                }

                // Fetch the sequence from the reference genome
                let seq = if is_multi {
                    get_seq_multi(&range, opt_genome)?
                } else {
                    get_seq(&range, opt_genome)?
                };

                //----------------------------
                // Output
                //----------------------------
                writer.write_all(format!(">{}\n{}\n", range, seq).as_ref())?;
            }

            // End of a block
            writer.write_all("\n".as_ref())?;
        }
    }

    Ok(())
}

fn get_seq(range: &intspan::Range, genome: &str) -> anyhow::Result<String> {
    let pos = format!("{}:{}-{}", range.chr(), range.start(), range.end());
    let mut gseq = intspan::get_seq_faidx(genome, &pos)?;

    if range.strand() == "-" {
        gseq = std::str::from_utf8(&bio::alphabets::dna::revcomp(gseq.bytes()))
            .unwrap()
            .to_string();
    }

    Ok(gseq)
}

fn get_seq_multi(range: &intspan::Range, genome: &str) -> anyhow::Result<String> {
    let pos = format!(
        "{}.{}:{}-{}",
        range.name(),
        range.chr(),
        range.start(),
        range.end()
    );
    let mut gseq = intspan::get_seq_faidx(genome, &pos)?;

    if range.strand() == "-" {
        gseq = std::str::from_utf8(&bio::alphabets::dna::revcomp(gseq.bytes()))
            .unwrap()
            .to_string();
    }

    Ok(gseq)
}
