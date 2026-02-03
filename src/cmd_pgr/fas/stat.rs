use clap::*;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("stat")
        .about("Basic statistics of block FA files")
        .after_help(
            r###"
Calculates basic statistics of block FA files (length, comparable, difference, etc.).

Input files can be gzipped. If the input file is 'stdin', data is read from standard input.

Examples:
1. Get statistics for block FA files:
   pgr fas stat tests/fas/part1.fas

2. Statistics treating the last sequence as an outgroup:
   pgr fas stat tests/fas/part1.fas --outgroup

3. Output results to a file:
   pgr fas stat tests/fas/part1.fas -o output.tsv

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
            Arg::new("has_outgroup")
                .long("outgroup")
                .action(ArgAction::SetTrue)
                .help("Indicates the presence of outgroups at the end of each block"),
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
    let has_outgroup = args.get_flag("has_outgroup");

    let field_names = vec![
        "target",
        "length",
        "comparable",
        "difference",
        "gap",
        "ambiguous",
        "D",
        "indel",
    ];

    //----------------------------
    // Operating
    //----------------------------
    writer.write_all(format!("{}\n", field_names.join("\t")).as_ref())?;

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = intspan::reader(infile);

        while let Ok(block) = pgr::libs::fas::next_fas_block(&mut reader) {
            let target = block.entries.first().unwrap().range().to_string();

            let mut seqs: Vec<&[u8]> = vec![];
            for entry in &block.entries {
                seqs.push(entry.seq().as_ref());
            }

            if has_outgroup {
                seqs.pop();
            }

            // let (length, comparable, difference, gap, ambiguous, mean_d) = alignment_stat(&seqs);
            let result = pgr::libs::alignment::alignment_stat(&seqs);

            let mut indel_ints = intspan::IntSpan::new();
            for seq in seqs {
                indel_ints.merge(&pgr::libs::alignment::indel_intspan(seq));
            }

            writer.write_all(
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    target,
                    result.0,
                    result.1,
                    result.2,
                    result.3,
                    result.4,
                    result.5,
                    indel_ints.span_size(),
                )
                .as_ref(),
            )?;
        }
    }

    Ok(())
}
