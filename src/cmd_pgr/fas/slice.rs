use clap::*;
use std::collections::BTreeMap;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("slice")
        .about("Extracts alignment slices")
        .after_help(
            r###"
Extracts alignment slices from block FA files using a runlist JSON.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* The JSON file (--required) keys are chromosome/sequence names, and values are runlists (e.g., "1-100,200-300")

Examples:
1. Extract slices defined in a JSON file:
   pgr fas slice tests/fas/slice.fas -r tests/fas/slice.json

2. Extract slices and name the output based on a specific species:
   pgr fas slice tests/fas/slice.fas -r tests/fas/slice.json --name S288c

3. Output results to a file:
   pgr fas slice tests/fas/slice.fas -r tests/fas/slice.json -o output.fas

"###,
        )
        .arg(
            Arg::new("runlist.json")
                .short('r')
                .long("required")
                .required(true)
                .num_args(1)
                .help("Required: JSON file describing ranges to extract"),
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
                .help("Reference species name. Default is the first species"),
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

    let json = intspan::read_json(args.get_one::<String>("runlist.json").unwrap());
    let set = intspan::json2set(&json);

    let mut name = if args.contains_id("name") {
        args.get_one::<String>("name").unwrap().to_string()
    } else {
        "".to_string()
    };

    //----------------------------
    // Operating
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = intspan::reader(infile);

        while let Ok(block) = pgr::libs::fas::next_fas_block(&mut reader) {
            // the first name of the first block
            if name.is_empty() {
                name = block.names.first().unwrap().to_string();
            }

            let idx = block.names.iter().position(|x| x == &name);
            if idx.is_none() {
                continue;
            }
            let trange = block.entries.get(idx.unwrap()).unwrap().range().clone();

            // chr present
            if !set.contains_key(trange.chr()) {
                continue;
            }
            if set.get(trange.chr()).unwrap().is_empty() {
                continue;
            }

            // has intersect
            let i_ints_chr = trange.intspan().intersect(set.get(trange.chr()).unwrap());
            if i_ints_chr.is_empty() {
                continue;
            }

            // target sequence intspan
            let t_ints_seq = pgr::libs::alignment::seq_intspan(
                block.entries.get(idx.unwrap()).unwrap().seq(),
            );

            // every sequence intspans
            let mut ints_seq_of = BTreeMap::new();
            // all indel region
            let mut indel_ints = intspan::IntSpan::new();
            for (i, name) in block.names.iter().enumerate() {
                let seq = block.entries.get(i).unwrap().seq();
                ints_seq_of.insert(
                    name.to_string(),
                    pgr::libs::alignment::seq_intspan(seq),
                );
                indel_ints.merge(&pgr::libs::alignment::indel_intspan(seq));
            }

            // there may be more than one subslice intersect this alignment
            let mut sub_slices: Vec<_> = vec![];
            for (lower, upper) in i_ints_chr.spans() {
                // chr positions to align
                let ss_start = pgr::libs::alignment::chr_to_align(
                    &t_ints_seq,
                    lower,
                    trange.start,
                    trange.strand(),
                )
                .unwrap();
                let ss_end = pgr::libs::alignment::chr_to_align(
                    &t_ints_seq,
                    upper,
                    trange.start,
                    trange.strand(),
                )
                .unwrap();
                if ss_start >= ss_end {
                    continue;
                }
                let mut ss_ints = intspan::IntSpan::from_pair(ss_start, ss_end);

                // borders of subslice inside an indel
                for n in [ss_start, ss_end] {
                    if indel_ints.contains(n) {
                        let island = indel_ints.find_islands_n(n);
                        ss_ints.subtract(&island);
                    }
                }
                sub_slices.push(ss_ints);
            }

            // write headers and sequences
            for ss in &sub_slices {
                let ss_start = ss.min();
                let ss_end = ss.max();

                // align positions to chromosomes of difference species
                for (i, name) in block.names.iter().enumerate() {
                    let range = block.entries.get(i).unwrap().range();
                    let start = pgr::libs::alignment::align_to_chr(
                        ints_seq_of.get(name).unwrap(),
                        ss_start,
                        range.start,
                        range.strand(),
                    )
                    .unwrap();
                    let end = pgr::libs::alignment::align_to_chr(
                        ints_seq_of.get(name).unwrap(),
                        ss_end,
                        range.start,
                        range.strand(),
                    )
                    .unwrap();
                    let ss_range = intspan::Range::from_full(
                        range.name(),
                        range.chr(),
                        range.strand(),
                        start,
                        end,
                    );

                    // seq of this sub slice
                    let ss_seq = &block.entries.get(i).unwrap().seq()
                        [(ss_start - 1) as usize..ss_end as usize];

                    //----------------------------
                    // Output
                    //----------------------------
                    writer.write_all(
                        format!(">{}\n{}\n", ss_range, std::str::from_utf8(ss_seq).unwrap())
                            .as_ref(),
                    )?;
                }
            }

            // end of a block
            writer.write_all("\n".as_ref())?;
        } // block
    } // infile

    Ok(())
}
