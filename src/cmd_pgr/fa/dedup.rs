use clap::*;
use std::collections::HashMap;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("dedup")
        .about("Deduplicates records in FASTA file(s)")
        .after_help(
            r###"
This command removes duplicate records from FASTA files.

Deduplication modes:
* By name (default): Compare sequence names only
* By description (-d): Compare full headers (name + description)
* By sequence (-s): Compare sequence contents

Comparison options:
* -b: Compare both strands (forward and reverse complement)
* -c: Case-insensitive comparison

Output options:
* -f FILE: Save duplicated entries mapping to FILE
* Format: original_name    duplicate_name

Notes:
* First occurrence is kept, others removed
* Supports both plain text and gzipped (.gz) files
* -b implies case-insensitive comparison for sequences

 sequence name
 | |
>sq0 LN:13
     |   |
     description

Examples:
1. Basic deduplication by name:
   pgr fa dedup input.fa -o output.fa

2. By sequence content:
   pgr fa dedup input.fa -s -o output.fa

3. Compare both strands:
   pgr fa dedup input.fa -s -b -o output.fa

4. Save duplicates mapping:
   pgr fa dedup input.fa -f dups.tsv -o output.fa

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input FASTA file(s) to process"),
        )
        .arg(
            Arg::new("desc")
                .long("desc")
                .short('d')
                .action(ArgAction::SetTrue)
                .help("Deduplicate by name and description"),
        )
        .arg(
            Arg::new("seq")
                .long("seq")
                .short('s')
                .action(ArgAction::SetTrue)
                .help("Deduplicate by sequence"),
        )
        .arg(
            Arg::new("both")
                .long("both")
                .short('b')
                .action(ArgAction::SetTrue)
                .help("Compare both strands (implies --case)"),
        )
        .arg(
            Arg::new("case")
                .long("case")
                .short('c')
                .action(ArgAction::SetTrue)
                .help("Case insensitive comparison"),
        )
        .arg(
            Arg::new("file")
                .long("file")
                .short('f')
                .num_args(1)
                .help("File to save duplicated names"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let is_desc = args.get_flag("desc");
    let is_seq = args.get_flag("seq");
    let is_both = args.get_flag("both");
    let is_insensitive = args.get_flag("case");

    let opts = pgr::libs::fasta::dedup::DedupOptions {
        is_seq,
        is_desc,
        is_both,
        is_insensitive,
    };

    let mut fa_out = pgr::libs::fmt::fa::writer(crate::cmd_pgr::args::get_outfile(args))?;

    //----------------------------
    // Process
    //----------------------------
    let mut subject_map: HashMap<u64, Vec<String>> = HashMap::new();

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut fa_in = pgr::libs::fmt::fa::reader(infile)?;

        for result in fa_in.records() {
            // obtain record or fail with error
            let record = result?;

            let name = record.name();
            let desc = record.description();
            let seq = record.sequence();

            let name_str = String::from_utf8(record.name().into())?;

            // filters
            let mut flag_pass = true;

            // name/desc/sequence to u64 signatures
            let subject = pgr::libs::fasta::dedup::record_signature(
                name,
                desc.map(|v| &**v),
                seq.get(..).unwrap(),
                &opts,
            )?;

            if let std::collections::hash_map::Entry::Vacant(e) = subject_map.entry(subject) {
                e.insert(vec![name_str]);
            } else {
                flag_pass = false;
                subject_map.get_mut(&subject).unwrap().push(name_str);
            }

            if !flag_pass {
                continue;
            }
            fa_out.write_record(&record)?;
        }
    }

    if args.contains_id("file") {
        let opt_file = args.get_one::<String>("file").unwrap();
        let mut writer = pgr::writer(opt_file)?;

        for v in subject_map.values() {
            if v.len() < 2 {
                continue;
            }

            for i in 1..v.len() {
                writer.write_fmt(format_args!("{}\t{}\n", v[0], v[i]))?;
            }
        }
    }

    Ok(())
}
