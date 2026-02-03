use clap::*;
use std::collections::BTreeMap;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("join")
        .about("Join multiple block fasta files by a common target")
        .after_help(
            r###"
Joins multiple block FA files based on a common target sequence.

Input files can be gzipped. If the input file is 'stdin', data is read from standard input.

Examples:
1. Join multiple block FA files:
   pgr fas join tests/fas/part1.fas tests/fas/part2.fas

2. Join files based on a specific species:
   pgr fas join tests/fas/part1.fas tests/fas/part2.fas --name S288c

3. Output results to a file:
   pgr fas join tests/fas/part1.fas tests/fas/part2.fas -o output.fas

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input block FA file(s) to process"),
        )
        .arg(
            Arg::new("name")
                .long("name")
                .num_args(1)
                .help("Target species name. Default is the first species"),
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

    let mut name = if args.contains_id("name") {
        args.get_one::<String>("name").unwrap().to_string()
    } else {
        "".to_string()
    };
    let mut block_of: BTreeMap<String, Vec<pgr::libs::fas::FasEntry>> = BTreeMap::new();

    //----------------------------
    // Operating
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = intspan::reader(infile);

        while let Ok(block) = pgr::libs::fas::next_fas_block(&mut reader) {
            if name.is_empty() {
                name = block.names.first().unwrap().to_string();
            }

            let idx = block.names.iter().position(|x| x == &name);
            if idx.is_none() {
                continue;
            }

            let idx = idx.unwrap();
            let header = block.entries.get(idx).unwrap().range().to_string();

            if !block_of.contains_key(&header) {
                // init
                block_of.insert(header.to_string(), vec![]);

                // entry with the selected name goes first
                block_of
                    .get_mut(&header)
                    .unwrap()
                    .push(block.entries.get(idx).unwrap().clone());
            }

            for entry in &block.entries {
                if entry.range().name() == &name {
                    continue;
                }
                block_of.get_mut(&header).unwrap().push(entry.clone());
            }
        }
    }

    //----------------------------
    // Output
    //----------------------------
    for v in block_of.values() {
        for e in v {
            writer.write_all(e.to_string().as_ref())?;
        }
        // end of a block
        writer.write_all("\n".as_ref())?;
    }

    Ok(())
}
