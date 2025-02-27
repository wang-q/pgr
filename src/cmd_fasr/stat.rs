use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("stat")
        .about("Extract a subset of species")
        .after_help(
            r###"
* <infiles> are paths to block fasta files, .fas.gz is supported
    * infile == stdin means reading from STDIN

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Set the input files to use"),
        )
        .arg(
            Arg::new("has_outgroup")
                .long("outgroup")
                .action(ArgAction::SetTrue)
                .help("There are outgroups at the end of each block"),
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

        while let Ok(block) = pgr::next_fas_block(&mut reader) {
            let target = block.entries.first().unwrap().range().to_string();

            let mut seqs: Vec<&[u8]> = vec![];
            for entry in &block.entries {
                seqs.push(entry.seq().as_ref());
            }

            if has_outgroup {
                seqs.pop();
            }

            // let (length, comparable, difference, gap, ambiguous, mean_d) = alignment_stat(&seqs);
            let result = pgr::alignment_stat(&seqs);

            let mut indel_ints = intspan::IntSpan::new();
            for seq in seqs {
                indel_ints.merge(&pgr::indel_intspan(seq));
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
