use anyhow::Context;
use clap::{Arg, ArgMatches, Command};
use pgr::libs::olc::layout::build_layouts;
use pgr::libs::olc::overlap::{Overlap, OverlapType};
use pgr::libs::paf::parser::parse_paf;
use std::collections::HashMap;
use std::io::{BufReader, Write};

/// Build the clap subcommand for layout.
pub fn make_subcommand() -> Command {
    Command::new("layout")
        .about("Chains unitigs into layouts from an overlap PAF (OLC stage 2)")
        .after_help(
            r###"
Builds greedy layouts from the exact overlaps produced by `pgr asm ovlp`:
every unitig end gets its best extension edge, unplaced unitigs are seeded
longest-first, and chains grow in both directions through mutual-best
junctions. Ambiguous junctions (two near-equal best partners, e.g. repeats)
and non-reciprocal edges stop the chain, so branches stay separate and no
heuristic picks a bubble path.

The unitig FASTA files must be the same files passed to `pgr asm ovlp` (the
`stem:name` prefixes are re-derived here and must match the PAF names).
The PAF file is the first positional argument.

Output is a layout TSV (no header), one line per step:
`contig_id<TAB>step<TAB>unitig_name<TAB>strand<TAB>q_start<TAB>q_end<TAB>overlap_len`
where q_start/q_end is the unitig's interval in the contig and overlap_len
is the exact overlap with the previous step (0 for the first step).

Examples:
1. Layout overlaps from two k values:
   pgr asm layout ovlp.paf k21.fa k51.fa -o layout.tsv
"###,
        )
        .arg(
            Arg::new("paf")
                .num_args(1)
                .index(1)
                .required(true)
                .help("Overlap PAF file (from pgr asm ovlp)"),
        )
        .arg(
            Arg::new("infiles")
                .num_args(1..)
                .index(2)
                .required(true)
                .help("Unitig FASTA file(s), same as passed to ovlp"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the layout command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let paf_path = args.get_one::<String>("paf").unwrap();
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    // Reject `-o` that would overwrite an input file (PAF or unitig FASTA).
    crate::cmd_pgr::args::ensure_outfile_distinct(
        outfile,
        infiles
            .iter()
            .map(|s| s.as_str())
            .chain(std::iter::once(paf_path.as_str())),
    )?;

    let unitigs = super::common::read_unitigs(&infiles)?;
    let mut id = HashMap::new();
    for (i, u) in unitigs.iter().enumerate() {
        id.insert(u.name.as_str(), i);
    }

    let reader = pgr::libs::io::reader(paf_path)
        .with_context(|| format!("failed to open input {paf_path}"))?;
    let records = parse_paf(BufReader::new(reader))?;
    let mut overlaps = Vec::with_capacity(records.len());
    for rec in &records {
        let qid = *id
            .get(rec.query_name.as_str())
            .with_context(|| format!("PAF query {} not found in unitigs", rec.query_name))?;
        let tid = *id
            .get(rec.target_name.as_str())
            .with_context(|| format!("PAF target {} not found in unitigs", rec.target_name))?;
        let otype = rec
            .tags
            .iter()
            .find_map(|t| t.strip_prefix("ov:A:"))
            .map(|c| {
                if c == "C" {
                    OverlapType::Contain
                } else {
                    OverlapType::Dovetail
                }
            })
            .unwrap_or(OverlapType::Dovetail);
        overlaps.push(Overlap {
            qid,
            tid,
            strand: rec.strand,
            q_start: rec.query_start as usize,
            q_end: rec.query_end as usize,
            t_start: rec.target_start as usize,
            t_end: rec.target_end as usize,
            length: rec.matches as usize,
            otype,
        });
    }

    let layouts = build_layouts(&unitigs, &overlaps)?;
    let mut out = pgr::libs::io::writer(outfile)
        .with_context(|| format!("failed to open output {outfile}"))?;
    super::common::write_layout_tsv(&mut out, &unitigs, &layouts)?;
    out.flush()?;
    Ok(())
}
