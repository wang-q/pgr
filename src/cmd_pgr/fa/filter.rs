use clap::*;
use std::collections::BTreeSet;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("filter")
        .about("Filters and formats sequences in FASTA file(s)")
        .after_help(
            r###"
This command filters and formats sequences in FASTA files.

Filters:
* --min-len N: Keep sequences >= N bp
* --max-len N: Keep sequences <= N bp
* --max-n N: Keep sequences with < N ambiguous bases
* --uniq: Remove duplicate sequence IDs

Formatters:
* --upper: Convert sequences to uppercase
* --iupac: Convert ambiguous codes to 'N'
* --dash: Remove dashes from sequences
* --simplify: Simplify sequence names (truncate at first space/./,/-)
* --line N: Set sequence line length

Notes:
        * Multiple filters can be combined
        * Supports both plain text and gzipped (.gz) files
        * For duplicate IDs, keeps the first occurrence
        * Not all faFilter options have been implemented
          Wildcards for names can be easily implemented with `pgr fa some`

        Examples:
        1. Filter by size:
           pgr fa filter input.fa --min-len 100 --max-len 1000

        2. Format sequences:
           pgr fa filter input.fa --upper --iupac --line 80

        3. Process multiple files:
           pgr fa filter *.fa --uniq --simplify -o output.fa

        "###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("FASTA"))
        .arg(
            Arg::new("min_len")
                .long("min-len")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Pass sequences at least this big"),
        )
        .arg(
            Arg::new("max_len")
                .long("max-len")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Pass sequences this size or smaller"),
        )
        .arg(
            Arg::new("max_n")
                .long("max-n")
                .short('n')
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Pass sequences with fewer than this number of Ns"),
        )
        .arg(
            Arg::new("uniq")
                .long("uniq")
                .short('u')
                .action(ArgAction::SetTrue)
                .help("Unique, removes duplicated ids, keeping the first"),
        )
        .arg(
            Arg::new("upper")
                .long("upper")
                .short('U')
                .action(ArgAction::SetTrue)
                .help("Convert all sequences to upper cases"),
        )
        .arg(
            Arg::new("iupac")
                .long("iupac")
                .short('N')
                .action(ArgAction::SetTrue)
                .help("Convert IUPAC ambiguous codes to 'N'"),
        )
        .arg(
            Arg::new("dash")
                .long("dash")
                .short('d')
                .action(ArgAction::SetTrue)
                .help("Remove dashes '-'"),
        )
        .arg(
            Arg::new("simplify")
                .long("simplify")
                .action(ArgAction::SetTrue)
                .help("Simplify sequence names"),
        )
        .arg(
            Arg::new("line")
                .long("line")
                .short('l')
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Sequence line length"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let opt_minsize = args
        .get_one::<usize>("min_len")
        .copied()
        .unwrap_or(pgr::libs::fasta::filter::NO_LIMIT);
    let opt_maxsize = args
        .get_one::<usize>("max_len")
        .copied()
        .unwrap_or(pgr::libs::fasta::filter::NO_LIMIT);
    let opt_maxn = args
        .get_one::<usize>("max_n")
        .copied()
        .unwrap_or(pgr::libs::fasta::filter::NO_LIMIT);
    let opt_line = args.get_one::<usize>("line").copied().unwrap_or(usize::MAX);

    let is_uniq = args.get_flag("uniq");
    let is_upper = args.get_flag("upper");
    let is_iupac = args.get_flag("iupac");
    let is_dash = args.get_flag("dash");
    let is_simplify = args.get_flag("simplify");

    let mut fa_out =
        pgr::libs::fmt::fa::writer_with_wrap(crate::cmd_pgr::args::get_outfile(args), opt_line)?;

    //----------------------------
    // Process
    //----------------------------
    let mut set_list: BTreeSet<String> = BTreeSet::new();
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut fa_in = pgr::libs::fmt::fa::reader(infile)?;

        for result in fa_in.records() {
            // obtain record or fail with error
            let record = result?;

            let mut name = String::from_utf8(record.name().into())?;
            if is_simplify {
                if let Some(i) = name.find(&[' ', '.', ',', '-'][..]) {
                    name = name[..i].to_string();
                }
            }
            let seq = record.sequence();

            // Apply filters
            if !pgr::libs::fasta::filter::pass_filters(
                seq.get(..).unwrap(),
                opt_minsize,
                opt_maxsize,
                opt_maxn,
                is_uniq,
                &mut set_list,
                &name,
            ) {
                continue;
            }

            // Apply formatters
            let seq_out = pgr::libs::fasta::filter::format_sequence(
                seq.get(..).unwrap(),
                is_dash,
                is_iupac,
                is_upper,
            );

            let record_out = pgr::libs::fmt::fa::new_record(&name, seq_out.as_bytes());
            fa_out.write_record(&record_out)?;
        }
    }

    Ok(())
}
