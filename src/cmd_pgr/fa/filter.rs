use anyhow::Context;
use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};
use std::collections::BTreeSet;
use std::io::Write;

/// Build the clap subcommand for filter.
pub fn make_subcommand() -> Command {
    Command::new("filter")
        .about("Filters and formats sequences in FASTA file(s)")
        .after_help(
            r###"
This command filters and formats sequences in FASTA files.

Filters:
* --min-len N: Keep sequences >= N bp
* --max-len N: Keep sequences <= N bp
* --max-n N: Keep sequences with <= N ambiguous bases (N/IUPAC)
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
* Reads from stdin if input file is 'stdin'
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
        .arg(crate::cmd_pgr::args::min_len_arg())
        .arg(crate::cmd_pgr::args::max_len_arg())
        .arg(
            Arg::new("max_n")
                .long("max-n")
                .short('n')
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Pass sequences with at most this number of ambiguous bases (N/IUPAC)"),
        )
        .arg(
            Arg::new("uniq")
                .long("uniq")
                .short('u')
                .action(ArgAction::SetTrue)
                .help("Unique, removes duplicated ids, keeping the first"),
        )
        .arg(crate::cmd_pgr::args::upper_arg())
        .arg(
            Arg::new("iupac")
                .long("iupac")
                .short('N')
                .action(ArgAction::SetTrue)
                .help("Convert IUPAC ambiguous codes to 'N'"),
        )
        .arg(crate::cmd_pgr::args::dash_arg())
        .arg(
            Arg::new("simplify")
                .long("simplify")
                .action(ArgAction::SetTrue)
                .help("Simplify sequence names"),
        )
        .arg(crate::cmd_pgr::args::line_arg(None))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the filter command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
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
    anyhow::ensure!(
        opt_line > 0,
        "--line must be positive (use a large value for no wrapping): {}",
        opt_line
    );

    let is_uniq = args.get_flag("uniq");
    let is_upper = args.get_flag("upper");
    let is_iupac = args.get_flag("iupac");
    let is_dash = args.get_flag("dash");
    let is_simplify = args.get_flag("simplify");

    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut fa_out = pgr::libs::fmt::fa::writer_with_wrap(outfile, opt_line)
        .with_context(|| format!("Failed to open writer for {}", outfile))?;

    let mut set_list: BTreeSet<String> = BTreeSet::new();
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut fa_in = pgr::libs::fmt::fa::reader(infile)
            .with_context(|| format!("Failed to open reader for {}", infile))?;

        for result in fa_in.records() {
            // obtain record or fail with error
            let record = result?;

            let mut name = String::from_utf8(record.name().into())?;
            if is_simplify {
                name = pgr::libs::io::simplify_name(&name).to_string();
            }
            let seq = record.sequence();

            // Apply filters
            if !pgr::libs::fasta::filter::pass_filters(
                seq.as_ref(),
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
                seq.as_ref(),
                is_dash,
                is_iupac,
                is_upper,
            );

            let record_out =
                pgr::libs::fmt::fa::new_record_preserving_desc(&name, &record, seq_out.as_bytes());
            fa_out.write_record(&record_out)?;
        }
    }

    fa_out.get_mut().flush()?;
    Ok(())
}
