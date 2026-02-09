use clap::*;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("consensus")
        .about("Generates consensus sequences using POA")
        .after_help(
            r###"
Generates consensus sequences using POA (Partial Order Alignment) graph.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* POA Engine:
    * `--engine builtin` (default): Uses built-in Rust implementation.
    * `--engine spoa`: Forces use of external `spoa` command.
* Alignment Parameters:
    * Configurable via `--match`, `--mismatch`, `--gap-open`, `--gap-extend`, `--algorithm`.
    * Defaults: Global alignment; Match 5, Mismatch -4, GapOpen -8, GapExtend -6.
* Supports parallel processing for improved performance
    * Running in parallel mode with 1 reader, 1 writer and the corresponding number of workers
    * The order of output may be different from the original
* If outgroups are present, they are handled appropriately

Examples:
1. Generate consensus sequences from a block FA file:
   pgr fas consensus tests/fas/example.fas

2. Generate consensus sequences using built-in engine:
   pgr fas consensus tests/fas/example.fas --engine builtin

3. Generate consensus sequences with outgroups:
   pgr fas consensus tests/fas/example.fas --outgroup

4. Run in parallel with 4 threads:
   pgr fas consensus tests/fas/example.fas --parallel 4

5. Output results to a file:
   pgr fas consensus tests/fas/example.fas -o output.fas

"###,
        )
        .arg(
            Arg::new("engine")
                .long("engine")
                .value_parser(["builtin", "spoa"])
                .default_value("builtin")
                .help("POA engine to use"),
        )
        .arg(
            Arg::new("match")
                .long("match")
                .short('m')
                .value_parser(value_parser!(i32))
                .default_value("5")
                .allow_negative_numbers(true)
                .help("Score for matching bases"),
        )
        .arg(
            Arg::new("mismatch")
                .long("mismatch")
                .short('n')
                .value_parser(value_parser!(i32))
                .default_value("-4")
                .allow_negative_numbers(true)
                .help("Score for mismatching bases"),
        )
        .arg(
            Arg::new("gap_open")
                .long("gap-open")
                .short('g')
                .value_parser(value_parser!(i32))
                .default_value("-8")
                .allow_negative_numbers(true)
                .help("Gap opening penalty"),
        )
        .arg(
            Arg::new("gap_extend")
                .long("gap-extend")
                .short('e')
                .value_parser(value_parser!(i32))
                .default_value("-6")
                .allow_negative_numbers(true)
                .help("Gap extension penalty"),
        )
        .arg(
            Arg::new("algorithm")
                .long("algorithm")
                .short('l')
                .value_parser(["local", "global", "semi_global"])
                .default_value("global") // Default to global for fas consensus
                .help("Alignment mode"),
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input block FA file(s) to process"),
        )
        .arg(
            Arg::new("cname")
                .long("cname")
                .num_args(1)
                .default_value("consensus")
                .help("Name of the consensus"),
        )
        .arg(
            Arg::new("has_outgroup")
                .long("outgroup")
                .action(ArgAction::SetTrue)
                .help("Indicates the presence of outgroups at the end of each block"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("1")
                .help("Number of threads for parallel processing"),
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
    let opt_parallel = *args.get_one::<usize>("parallel").unwrap();

    //----------------------------
    // Operating
    //----------------------------
    if opt_parallel == 1 {
        // Single-threaded mode
        let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap());

        for infile in args.get_many::<String>("infiles").unwrap() {
            let mut reader = pgr::reader(infile);
            while let Ok(block) = pgr::libs::fas::next_fas_block(&mut reader) {
                let out_string = proc_block(&block, args)?;
                writer.write_all(out_string.as_ref())?;
            }
        }
    } else {
        // Parallel mode
        proc_block_p(args)?;
    }

    Ok(())
}

fn proc_block(block: &pgr::libs::fas::FasBlock, args: &ArgMatches) -> anyhow::Result<String> {
    //----------------------------
    // Args
    //----------------------------
    let cname = args.get_one::<String>("cname").unwrap();
    let has_outgroup = args.get_flag("has_outgroup");

    let engine = args.get_one::<String>("engine").unwrap();

    let match_score = *args.get_one::<i32>("match").unwrap();
    let mismatch_score = *args.get_one::<i32>("mismatch").unwrap();
    let gap_open = *args.get_one::<i32>("gap_open").unwrap();
    let gap_extend = *args.get_one::<i32>("gap_extend").unwrap();
    let algorithm = args.get_one::<String>("algorithm").unwrap();

    // Map algorithm string to integer code (0=local, 1=global, 2=semi_global) for internal use/spoa
    let algo_code = match algorithm.as_str() {
        "local" => 0,
        "global" => 1,
        "semi_global" => 2,
        _ => 1,
    };

    //----------------------------
    // Ops
    //----------------------------
    let mut seqs = vec![];

    let outgroup = if has_outgroup {
        Some(block.entries.iter().last().unwrap())
    } else {
        None
    };

    for entry in &block.entries {
        seqs.push(entry.seq().as_ref());
    }
    if outgroup.is_some() {
        seqs.pop().unwrap(); // Remove the outgroup sequence
    }

    // Generate consensus sequence
    let mut cons = match engine.as_str() {
        "spoa" => pgr::libs::alignment::get_consensus_poa_external(
            &seqs,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
            algo_code,
        )
        .unwrap(),
        "builtin" | _ => pgr::libs::alignment::get_consensus_poa_builtin(
            &seqs,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
            algo_code,
        )
        .unwrap(),
    };
    cons = cons.replace('-', "");

    let mut range = block.entries.first().unwrap().range().clone();

    //----------------------------
    // Output
    //----------------------------
    let mut out_string = "".to_string();
    if range.is_valid() {
        *range.name_mut() = cname.to_string();
        out_string += format!(">{}\n{}\n", range, cons).as_ref();
    } else {
        out_string += format!(">{}\n{}\n", cname, cons).as_ref();
    }
    if outgroup.is_some() {
        out_string += outgroup.unwrap().to_string().as_ref();
    }

    // end of a block
    out_string += "\n";

    Ok(out_string)
}

// Adopt from https://rust-lang-nursery.github.io/rust-cookbook/concurrency/threads.html#create-a-parallel-pipeline
fn proc_block_p(args: &ArgMatches) -> anyhow::Result<()> {
    let parallel = *args.get_one::<usize>("parallel").unwrap();
    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap());

    // Channel 1 - Read files to blocks
    let (snd1, rcv1) = crossbeam::channel::bounded::<pgr::libs::fas::FasBlock>(10);
    // Channel 2 - Results
    let (snd2, rcv2) = crossbeam::channel::bounded(10);

    crossbeam::scope(|s| {
        //----------------------------
        // Reader thread
        //----------------------------
        s.spawn(|_| {
            for infile in args.get_many::<String>("infiles").unwrap() {
                let mut reader = pgr::reader(infile);
                while let Ok(block) = pgr::libs::fas::next_fas_block(&mut reader) {
                    snd1.send(block).unwrap();
                }
            }
            // Close the channel - this is necessary to exit the for-loop in the worker
            drop(snd1);
        });

        //----------------------------
        // Worker threads
        //----------------------------
        for _ in 0..parallel {
            // Send to sink, receive from source
            let (sendr, recvr) = (snd2.clone(), rcv1.clone());
            // Spawn workers in separate threads
            s.spawn(move |_| {
                // Receive until channel closes
                for block in recvr.iter() {
                    let out_string = proc_block(&block, args).unwrap();
                    sendr.send(out_string).unwrap();
                }
            });
        }
        // Close the channel, otherwise sink will never exit the for-loop
        drop(snd2);

        //----------------------------
        // Writer thread
        //----------------------------
        for out_string in rcv2.iter() {
            writer.write_all(out_string.as_ref()).unwrap();
        }
    })
    .unwrap();

    Ok(())
}
