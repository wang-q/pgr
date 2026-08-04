//! `pgr rg span` — line-level span operations on `.rg` lines.

use clap::{builder, value_parser, Arg, ArgAction, ArgMatches, Command};
use pgr::libs::ds::Range;
use std::io::{BufRead, Write};

/// Clamp a valid result back into the representable `.rg` coordinate domain
/// (`1..=POS_INF - 1`). Saturating arithmetic in the Range ops can leave a
/// valid-looking range above the maximum (e.g. `shift` by a huge `-n`); such
/// output could never be re-read by the other rg commands, so it is clamped
/// instead of silently dropped downstream.
fn clamp_to_domain(range: Range) -> Range {
    if !range.is_valid() {
        return range;
    }
    let max = pgr::libs::ds::IntSpan::new().get_pos_inf();
    let start = *range.start().min(&max).max(&1);
    let end = *range.end().min(&max);
    if start <= end {
        Range::from_full(range.name(), range.chr(), range.strand(), start, end)
    } else {
        // Collapse to the invalid (0, 0) form, keeping the line identity.
        Range::from_full(range.name(), range.chr(), range.strand(), 0, 0)
    }
}

/// Build the clap subcommand for span.
pub fn make_subcommand() -> Command {
    Command::new("span")
        .about("Operates on spans in .rg files")
        .after_help(
            r###"
Applies an operation to each `.rg` line and writes the new range (or appends
it with `--append`). Operations:

* trim: remove N integers from the ends (5p / 3p / both)
* pad: add N integers to the ends
* shift: shift the range N bases toward the 5p or 3p end
* flank: retrieve the flanking region of size N at the 5p or 3p end
* excise: drop ranges smaller than N (written as an empty line)

Examples:
1. Trim 10 bp from both ends:
   pgr rg span a.rg --op trim -n 10
2. Flank -1 bp at the 3p end, appending the new range:
   pgr rg span a.rg --op flank -m 3p -n=-1 -a
"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input .rg files to process"),
        )
        .arg(
            Arg::new("op")
                .long("op")
                .num_args(1)
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("trim"),
                    builder::PossibleValue::new("pad"),
                    builder::PossibleValue::new("shift"),
                    builder::PossibleValue::new("flank"),
                    builder::PossibleValue::new("excise"),
                ])
                .default_value("trim")
                .help("Operation to perform"),
        )
        .arg(
            Arg::new("mode")
                .long("mode")
                .short('m')
                .num_args(1)
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("both"),
                    builder::PossibleValue::new("5p"),
                    builder::PossibleValue::new("3p"),
                ])
                .default_value("both")
                .help("Mode of the operation (5p/3p required for shift/flank)"),
        )
        .arg(
            Arg::new("number")
                .long("number")
                .short('n')
                .num_args(1)
                .value_parser(value_parser!(i32))
                .default_value("0")
                .help("Number of bases to trim, pad, shift or flank; length threshold for excise"),
        )
        .arg(
            Arg::new("append")
                .long("append")
                .short('a')
                .action(ArgAction::SetTrue)
                .help("Append the new range to the original line"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the span command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let op = args.get_one::<String>("op").unwrap().as_str();
    let mode = args.get_one::<String>("mode").unwrap().as_str();
    let number = *args.get_one::<i32>("number").unwrap();
    let is_append = args.get_flag("append");

    // Validate the op/mode combination before reading any input, so an
    // invalid invocation fails even when the input files are empty.
    if matches!(op, "shift" | "flank") && mode == "both" {
        anyhow::bail!("--mode both is invalid for {op}");
    }
    crate::cmd_pgr::args::ensure_outfile_distinct(
        outfile,
        args.get_many::<String>("infiles")
            .unwrap()
            .map(String::as_str),
    )?;

    let mut writer = pgr::writer(outfile)?;
    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = pgr::reader(infile)?;
        for line in reader.lines() {
            let line = line?;
            if line.trim_start().starts_with('#') {
                continue;
            }
            let range = pgr::libs::ds::Range::from_str(&line);
            if !pgr::libs::runlist::usable_range(&range) {
                continue;
            }
            let new: Range = match op {
                "trim" => match mode {
                    "5p" => range.trim_5p(number),
                    "3p" => range.trim_3p(number),
                    _ => range.trim(number),
                },
                "pad" => match mode {
                    // `saturating_neg` keeps `-n i32::MIN` from overflowing
                    // before the (already saturating) Range ops run.
                    "5p" => range.trim_5p(number.saturating_neg()),
                    "3p" => range.trim_3p(number.saturating_neg()),
                    _ => range.trim(number.saturating_neg()),
                },
                "shift" => match mode {
                    "5p" => range.shift_5p(number),
                    "3p" => range.shift_3p(number),
                    _ => unreachable!("mode validated before reading input"),
                },
                "flank" => match mode {
                    "5p" => range.flank_5p(number),
                    "3p" => range.flank_3p(number),
                    _ => unreachable!("mode validated before reading input"),
                },
                "excise" => {
                    if range.intspan().size() >= number {
                        range.clone()
                    } else {
                        Range::new()
                    }
                }
                _ => unreachable!("invalid span op"),
            };
            let new = clamp_to_domain(new);
            if is_append {
                writeln!(writer, "{}\t{}", line, new)?;
            } else {
                writeln!(writer, "{}", new)?;
            }
        }
    }
    writer.flush()?;
    Ok(())
}
