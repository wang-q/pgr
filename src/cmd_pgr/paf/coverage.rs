use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use pgr::libs::paf::cov::coverage_segments;
use pgr::libs::paf::parser::parse_paf;
use std::io::{BufReader, Write};

/// Build the clap subcommand for coverage.
pub fn make_subcommand() -> Command {
    Command::new("coverage")
        .about("Computes per-target coverage depth from PAF cg:Z tags")
        .after_help(
            r###"
Accumulates per-target alignment depth from the `cg:Z` CIGAR tag of every
PAF record (M/=/X/D operations cover the target, insertions do not) and
writes maximal segments of constant depth at or above `--minimum` as TSV:
`target<TAB>start<TAB>end<TAB>depth` (0-based half-open, sorted by target
then start).

Records without a `cg:Z` tag (e.g. `pgr psl to-paf` output) contribute
nothing. Pipelines that carry the tag include `pgr pl chainnet` →
`pgr maf to-paf`.

Examples:
1. All covered segments:
   pgr paf coverage ovlp.paf -o cov.tsv
2. Segments with depth at least 5:
   pgr paf coverage ovlp.paf -m 5 -o cov.tsv
"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg_with_numargs(
            "Input PAF file(s)",
            1..,
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("minimum")
                .long("minimum")
                .short('m')
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(u32))
                .help("Minimum depth to report"),
        )
}

/// Execute the coverage command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let minimum = *args.get_one::<u32>("minimum").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    let mut records = Vec::new();
    for infile in &infiles {
        let reader = pgr::libs::io::reader(infile)
            .with_context(|| format!("failed to open input {infile}"))?;
        records.extend(parse_paf(BufReader::new(reader))?);
    }
    let segs = coverage_segments(&records, minimum)?;
    let mut out = pgr::libs::io::writer(outfile)
        .with_context(|| format!("failed to open output {outfile}"))?;
    for (target, segs) in &segs {
        for (start, end, depth) in segs {
            writeln!(out, "{target}\t{start}\t{end}\t{depth}")?;
        }
    }
    out.flush()?;
    Ok(())
}
