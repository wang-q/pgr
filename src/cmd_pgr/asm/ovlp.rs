use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use pgr::libs::olc::overlap::{find_overlaps, OverlapOptions};
use std::io::Write;

/// Build the clap subcommand for ovlp.
pub fn make_subcommand() -> Command {
    Command::new("ovlp")
        .about("Finds exact overlaps between unitigs (OLC stage 1)")
        .after_help(
            r###"
Finds exact suffix/prefix overlaps between unitigs by seeding a canonical
k-mer index with the boundary k-mers of every unitig and verifying each
candidate by extension, so overlaps are exact and error-free (unitigs come
from the de Bruijn graph). This is the overlap stage of the OLC assembly
pipeline (see notes/design/olc.md); the caller is expected to assemble
unitigs at several k values first and pass the FASTA files here.

Overlaps are written as PAF with an `ov:A:D` (dovetail) or `ov:A:C`
(contain) tag. Unitig names are prefixed with the input file stem
(`stem:name`) so identical `unitig_<id>` names across k files stay unique;
the prefix is deterministic (only `[A-Za-z0-9_.-]` are kept).

Notes:
* Seed k is clamped to the shortest unitig length; unitigs shorter than the
  seed still appear as overlap targets
* Self overlaps, reverse-complement self matches, and overlaps below
  --min-overlap are discarded
* Output is sorted and deterministic

Examples:
1. Overlap unitigs from two k values:
   pgr asm ovlp k21.fa k51.fa -o ovlp.paf
2. Raise the seed and minimum overlap:
   pgr asm ovlp unitigs.fa -o ovlp.paf --overlap-k 21 --min-overlap 51
"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg_with_numargs(
            "Unitig FASTA file(s) to compare",
            1..,
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("overlap_k")
                .long("overlap-k")
                .num_args(1)
                .default_value("17")
                .value_parser(value_parser!(usize))
                .help("Seed k-mer length (clamped to the shortest unitig)"),
        )
        .arg(
            Arg::new("min_overlap")
                .long("min-overlap")
                .num_args(1)
                .default_value("34")
                .value_parser(value_parser!(usize))
                .help("Minimum accepted overlap length in bases"),
        )
}

/// Execute the ovlp command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let seed_k = *args.get_one::<usize>("overlap_k").unwrap();
    let min_overlap = *args.get_one::<usize>("min_overlap").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    // Reject `-o` that would overwrite an input file (unitig FASTA).
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, infiles.iter().map(|s| s.as_str()))?;

    let unitigs = super::common::read_unitigs(&infiles)?;
    let overlaps = find_overlaps(
        &unitigs,
        &OverlapOptions {
            seed_k,
            min_overlap,
        },
    )?;
    let mut out = pgr::libs::io::writer(outfile)
        .with_context(|| format!("failed to open output {outfile}"))?;
    for ov in &overlaps {
        let rec = super::common::to_paf(ov, &unitigs);
        pgr::libs::paf::record::write_paf_record(&mut out, &rec)?;
    }
    out.flush()?;
    Ok(())
}
