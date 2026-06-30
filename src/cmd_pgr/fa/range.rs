use clap::*;
use pgr::libs::loc;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("range")
        .about("Extracts sequence regions by coordinates")
        .after_help(
            r###"
This command extracts sequence regions from FASTA files using genomic coordinates.

Range format:
    seq_name(strand):start-end

Notes:
* Cannot read from stdin or plain gzip (requires BGZF for random access)
* Supports BGZF compressed files (.gz)
* Automatic index creation (.loc)
* LRU caching for better performance
* Reverse complement for negative strand
* All coordinates (<start> and <end>) are based on the positive strand
* Sort range file for better performance
* Cache size affects memory usage

Examples:
1. Single range:
   pgr fa range input.fa "chr1:1-1000"

2. Multiple ranges:
   pgr fa range input.fa "chr1:1-1000" "chr2(-):2000-3000"

3. From range file with larger cache:
   pgr fa range input.fa -r ranges.txt -c 10

4. Force update the index file:
   pgr fa range input.fa "chr1:1-1000" --update

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input FASTA file to process"),
        )
        .arg(
            Arg::new("ranges")
                .required(false)
                .index(2)
                .num_args(0..)
                .help("Ranges of interest"),
        )
        .arg(crate::cmd_pgr::args::rgfile_arg())
        .arg(
            Arg::new("cache")
                .long("cache")
                .short('c')
                .value_parser(value_parser!(std::num::NonZeroUsize))
                .num_args(1)
                .default_value("1")
                .help("Set the capacity of the LRU cache"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("update")
                .long("update")
                .short('u')
                .action(ArgAction::SetTrue)
                .help("Force update the .loc index file"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let infile = args.get_one::<String>("infile").unwrap();

    let mut fa_out = pgr::libs::fmt::fa::writer(args.get_one::<String>("outfile").unwrap())?;

    let ranges = crate::cmd_pgr::args::collect_ranges(args)?;

    let opt_cache = *args.get_one::<std::num::NonZeroUsize>("cache").unwrap();
    let mut cache: lru::LruCache<String, noodles_fasta::Record> = lru::LruCache::new(opt_cache);

    //----------------------------
    // Open files
    //----------------------------
    let force_update = args.get_flag("update");
    let (mut reader, loc_of) = loc::open_indexed(infile, force_update)?;

    //----------------------------
    // Output
    //----------------------------
    for el in ranges.iter() {
        let rg = intspan::Range::from_str(el);
        let seq_id = rg.chr().to_string();
        if !loc_of.contains_key(&seq_id) {
            log::warn!("{} for [{}] not found in the .loc index file", seq_id, el);
            continue;
        }

        if !cache.contains(&seq_id) {
            let record = loc::fetch_record(&mut reader, &loc_of, &seq_id)?;
            cache.put(seq_id.clone(), record);
        }

        let record: &noodles_fasta::Record = cache.get(&seq_id).unwrap();

        // name only
        if *rg.start() == 0 {
            fa_out.write_record(record)?;
            continue;
        }

        let definition = noodles_fasta::record::Definition::new(rg.to_string(), None);

        // slice here is 1-based
        let start = noodles_core::Position::new(*rg.start() as usize).unwrap();
        let end = noodles_core::Position::new(*rg.end() as usize).unwrap();

        let mut slice = record.sequence().slice(start..=end).unwrap();
        if rg.strand() == "-" {
            slice = slice.complement().rev().collect::<Result<_, _>>()?;
        }
        let record_rg = noodles_fasta::Record::new(definition, slice);

        fa_out.write_record(&record_rg)?;
    }

    Ok(())
}

// fn print_type_of<T: ?Sized>(_: &T) {
//     println!("{}", std::any::type_name::<T>())
// }
