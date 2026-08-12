use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use pgr::libs::olc::consensus::consensus;
use pgr::libs::olc::layout::{Layout, LayoutStep};
use std::collections::HashMap;
use std::io::{BufRead, Write};

/// Build the clap subcommand for cns.
pub fn make_subcommand() -> Command {
    Command::new("cns")
        .about("Stitches layouts into consensus contigs (OLC stage 3)")
        .after_help(
            r###"
Stitches the layouts produced by `pgr asm layout` into consensus contigs.
Overlaps are exact, so each layout is walked in order, every unitig is
oriented by its strand, and only the bases beyond the exact overlap with
the previous step are appended. A layout whose overlapping bases disagree
with the already-stitched contig is reported as an error (exact overlaps
must agree).

The unitig FASTA files must be the same files passed to `pgr asm ovlp` and
`pgr asm layout` (the `stem:name` prefixes are re-derived here and must
match). The layout TSV is the first positional argument.

Output is FASTA (`>contig_<id>,len=...,cov=...`, 70-column wrap, longest
first). `cov` is the approximate unitig depth (sum of unitig lengths over
the contig length).

Examples:
1. Consensus from a layout:
   pgr asm cns layout.tsv k21.fa k51.fa -o contigs.fa
2. Drop short contigs:
   pgr asm cns layout.tsv unitigs.fa -o contigs.fa --min-contig-len 500
"###,
        )
        .arg(
            Arg::new("layout")
                .num_args(1)
                .index(1)
                .required(true)
                .help("Layout TSV file (from pgr asm layout)"),
        )
        .arg(
            Arg::new("infiles")
                .num_args(1..)
                .index(2)
                .required(true)
                .help("Unitig FASTA file(s), same as passed to ovlp"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("min_contig_len")
                .long("min-contig-len")
                .num_args(1)
                .default_value("500")
                .value_parser(value_parser!(usize))
                .help("Minimum contig length in bases"),
        )
}

/// Execute the cns command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let layout_path = args.get_one::<String>("layout").unwrap();
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let min_contig_len = *args.get_one::<usize>("min_contig_len").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    // Reject `-o` that would overwrite an input file (layout TSV or unitig
    // FASTA).
    crate::cmd_pgr::args::ensure_outfile_distinct(
        outfile,
        infiles
            .iter()
            .map(|s| s.as_str())
            .chain(std::iter::once(layout_path.as_str())),
    )?;

    let unitigs = super::common::read_unitigs(&infiles)?;
    let mut id = HashMap::new();
    for (i, u) in unitigs.iter().enumerate() {
        id.insert(u.name.as_str(), i);
    }
    let layouts = parse_layouts(layout_path, &id)?;
    let contigs = consensus(&unitigs, &layouts, min_contig_len)?;

    let mut out = pgr::libs::io::writer(outfile)
        .with_context(|| format!("failed to open output {outfile}"))?;
    for (i, c) in contigs.iter().enumerate() {
        writeln!(
            out,
            ">contig_{},len={},cov={}",
            i + 1,
            c.seq.len(),
            super::common::format_cov(c.coverage)
        )?;
        for chunk in c.seq.chunks(70) {
            out.write_all(chunk)?;
            out.write_all(b"\n")?;
        }
    }
    out.flush()?;
    Ok(())
}

/// Parses a layout TSV into layouts with unitig ids.
fn parse_layouts(path: &str, id: &HashMap<&str, usize>) -> anyhow::Result<Vec<Layout>> {
    let reader =
        pgr::libs::io::reader(path).with_context(|| format!("failed to open input {path}"))?;
    let mut layouts: Vec<Layout> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        anyhow::ensure!(
            fields.len() == 7,
            "invalid layout line, expected 7 fields: {line}"
        );
        // `contig_N` ids are 1-based; reject `contig_0` (would underflow).
        let ci = {
            let n = parse_index(fields[0], "contig id")?;
            anyhow::ensure!(n >= 1, "invalid contig id in layout line: {line}");
            n - 1
        };
        let si = parse_index(fields[1], "step")?;
        let name = fields[2];
        let strand = match fields[3] {
            "+" => '+',
            "-" => '-',
            _ => anyhow::bail!("invalid strand in layout line: {line}"),
        };
        let q_start = fields[4]
            .parse::<usize>()
            .with_context(|| format!("invalid q_start in layout line: {line}"))?;
        let q_end = fields[5]
            .parse::<usize>()
            .with_context(|| format!("invalid q_end in layout line: {line}"))?;
        let overlap_len = fields[6]
            .parse::<usize>()
            .with_context(|| format!("invalid overlap_len in layout line: {line}"))?;
        let unitig = *id
            .get(name)
            .with_context(|| format!("layout unitig {name} not found in unitigs"))?;
        if layouts.len() <= ci {
            layouts.resize(ci + 1, Layout { steps: Vec::new() });
        }
        anyhow::ensure!(
            layouts[ci].steps.len() == si,
            "layout steps must be contiguous, got step {si} for {name}"
        );
        layouts[ci].steps.push(LayoutStep {
            unitig,
            strand,
            q_start,
            q_end,
            overlap_len,
        });
    }
    Ok(layouts)
}

/// Parses `contig_N` / step integers.
fn parse_index(s: &str, what: &str) -> anyhow::Result<usize> {
    let n = s.rsplit_once('_').map(|(_, n)| n).unwrap_or(s);
    n.parse::<usize>()
        .with_context(|| format!("invalid {what}: {s}"))
}
