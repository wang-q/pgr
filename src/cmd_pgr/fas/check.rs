use clap::*;
use indexmap::IndexMap;
use pgr::libs::loc;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("check")
        .about("Checks genome locations in block FA headers")
        .after_help(
            r###"
Checks genome locations in block FA headers against a chrom.sizes file.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

"###,
        )
        .arg(
            Arg::new("genome.fa")
                .short('r')
                .long("genome")
                .required(true)
                .num_args(1)
                .help("Path to the reference genome FA file"),
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input block FA file(s) to check"),
        )
        .arg(
            Arg::new("name")
                .long("name")
                .num_args(1)
                .help("Check sequences for a specific species"),
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

    //----------------------------
    // Ops
    //----------------------------
    let is_bgzf = {
        let path = std::path::Path::new(opt_genome);
        path.extension() == Some(std::ffi::OsStr::new("gz"))
    };
    let loc_file = format!("{}.loc", opt_genome);
    if !std::path::Path::new(&loc_file).is_file() {
        loc::create_loc(opt_genome, &loc_file, is_bgzf)?;
    }
    let loc_of: IndexMap<String, (u64, usize)> = loc::load_loc(&loc_file)?;

    let mut genome_reader = if is_bgzf {
        loc::Input::Bgzf(
            noodles_bgzf::io::indexed_reader::Builder::default().build_from_path(opt_genome)?,
        )
    } else {
        loc::Input::File(std::fs::File::open(std::path::Path::new(opt_genome))?)
    };

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = intspan::reader(infile);

        while let Ok(block) = pgr::next_fas_block(&mut reader) {
            let block_names = block.names;

            // Check if a specific species is requested
            if !opt_name.is_empty() && block_names.contains(opt_name) {
                for entry in &block.entries {
                    let entry_name = entry.range().name();
                    if entry_name == opt_name {
                        let status = check_seq(entry, &mut genome_reader, &loc_of)?;
                        writer.write_all(format!("{}\t{}\n", entry.range(), status).as_ref())?;
                    }
                }
            } else if opt_name.is_empty() {
                // Check all sequences in the block
                for entry in &block.entries {
                    let status = check_seq(entry, &mut genome_reader, &loc_of)?;
                    writer.write_all(format!("{}\t{}\n", entry.range(), status).as_ref())?;
                }
            }
        }
    }

    Ok(())
}

fn check_seq(
    entry: &pgr::FasEntry,
    reader: &mut loc::Input,
    loc_of: &IndexMap<String, (u64, usize)>,
) -> anyhow::Result<String> {
    let range = entry.range();
    let seq = entry.seq().to_vec();
    let seq = std::str::from_utf8(&seq)?
        .to_string()
        .to_ascii_uppercase()
        .replace('-', "");

    let gseq = if loc_of.contains_key(range.chr()) {
        loc::fetch_range_seq(reader, loc_of, range)?.to_ascii_uppercase()
    } else {
        String::new()
    };

    let status = if seq == gseq { "OK" } else { "FAILED" };

    Ok(status.to_string())
}
