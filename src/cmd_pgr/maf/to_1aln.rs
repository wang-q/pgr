//! `pgr maf to-1aln`: convert a contig-level MAF file into FastGA `.1aln`.
//!
//! This is the write-side sibling of `pgr paf to-1aln`, following the design
//! doc `notes/design/1aln.md` §6.1.1: a two-sequence MAF block carries the
//! base-level path (`s` line texts with explicit gaps) and enough skeleton
//! information (`srcSize`) to build the `.1aln` GDB skeleton without source
//! genomes.
//!
//! Axis convention: the `.1aln` `a` side is always forward. The first MAF
//! component is used as `a` when it is on the `+` strand; otherwise the second
//! is used and the first becomes the `b` side. A `-` strand `b` side is marked
//! reverse-complemented (`R` line) and its coordinates are reversed onto the
//! forward source.

use anyhow::{anyhow, bail, Context, Result};
use clap::{ArgMatches, Command};
use std::collections::HashMap;
use std::io::BufReader;

use pgr::libs::alignment::coords::reverse_range_pair;
use pgr::libs::fmt::maf::{next_maf_block, MafComp};
use pgr::libs::onepack::write::{
    cigar2tp, open_aln_writer, write_aln_record, write_skeleton_contigs, write_tspace,
};
use pgr::libs::paf::cigar::{cigar_from_alignment, reverse_cigar};

/// Build the clap subcommand for to-1aln.
pub fn make_subcommand() -> Command {
    Command::new("to-1aln")
        .about("Converts a contig-level MAF file to FastGA .1aln format")
        .after_help(
            r###"
Reads a contig-level MAF (two-sequence blocks where each `s` line is a top-level
source contig) and writes a FastGA `.1aln` (ONEcode trace-point) file. The MAF
`s` line texts carry the base-level path (explicit gaps), and `srcSize` supplies
the skeleton, so no source genomes are needed.

The `.1aln` `a` side is always forward. The first MAF component is used as `a`
when it is on the `+` strand; otherwise the second component is used and the
first becomes the `b` side. A `-` strand `b` side is marked reverse-complemented
(`R` line) with its coordinates reversed onto the forward source.

Multi-sequence MAF blocks (more than two `s` lines) are skipped with a warning.

Notes:
* The ONEcode container is binary and requires a real output file (not stdout).
* The skeleton is built from the MAF `s` line `src` names and `srcSize`; no
  source FASTA is required.

Examples:
1. Convert a contig-level MAF to .1aln:
   pgr maf to-1aln align.maf -o align.1aln
2. Use a custom trace-point spacing:
   pgr maf to-1aln align.maf -o align.1aln --tspace 50

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input MAF file",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg_required())
        .arg(
            clap::Arg::new("tspace")
                .long("tspace")
                .short('t')
                .num_args(1)
                .default_value("100")
                .value_parser(clap::value_parser!(i64))
                .help("Trace point spacing (default: 100)"),
        )
}

/// Execute the to-1aln command.
pub fn execute(args: &ArgMatches) -> Result<()> {
    let infile = args
        .get_one::<String>("infile")
        .context("missing required argument: infile")?;
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, [infile.as_str()])?;
    let tspace = *args.get_one::<i64>("tspace").unwrap();

    if outfile == "stdout" {
        bail!("to-1aln requires a real output file (the ONEcode container is binary)");
    }

    // First pass: collect the skeleton (contig names/lengths) from all blocks.
    let mut a_seqs: Vec<(String, i64)> = Vec::new();
    let mut b_seqs: Vec<(String, i64)> = Vec::new();
    let mut a_index: HashMap<String, usize> = HashMap::new();
    let mut b_index: HashMap<String, usize> = HashMap::new();
    {
        let file = std::fs::File::open(infile)
            .with_context(|| format!("Failed to open MAF file {infile}"))?;
        let mut reader = BufReader::new(file);
        while let Ok(block) = next_maf_block(&mut reader) {
            let Some((a, b)) = orient_pair(&block) else {
                continue;
            };
            insert_contig(&mut a_seqs, &mut a_index, a);
            insert_contig(&mut b_seqs, &mut b_index, b);
        }
    }
    if a_seqs.is_empty() {
        bail!("no usable two-sequence MAF blocks found in {infile}");
    }

    let mut writer = open_aln_writer(outfile)?;
    writer.add_reference(infile, 1);
    write_tspace(&mut writer, tspace)?;
    write_skeleton_contigs(&mut writer, &a_seqs)?;
    write_skeleton_contigs(&mut writer, &b_seqs)?;

    // Second pass: write one `.1aln` record per block.
    {
        let file = std::fs::File::open(infile)
            .with_context(|| format!("Failed to open MAF file {infile}"))?;
        let mut reader = BufReader::new(file);
        while let Ok(block) = next_maf_block(&mut reader) {
            let Some((a, b)) = orient_pair(&block) else {
                eprintln!("Warning: skipping non-two-sequence MAF block");
                continue;
            };
            write_block(&mut writer, a, b, &a_index, &b_index, tspace)?;
        }
    }

    writer.close()?;
    Ok(())
}

/// Orient the two components of a block so `a` is forward and `b` is the other.
///
/// Returns `None` for blocks that are not two-sequence (skipped by the caller).
fn orient_pair(block: &pgr::libs::fmt::maf::MafAli) -> Option<(&MafComp, &MafComp)> {
    if block.components.len() != 2 {
        return None;
    }
    let (c0, c1) = (&block.components[0], &block.components[1]);
    if c0.strand == '+' {
        Some((c0, c1))
    } else if c1.strand == '+' {
        Some((c1, c0))
    } else {
        // Both reverse: cannot orient a forward `a` side.
        None
    }
}

/// Insert a contig into the skeleton if not already present.
fn insert_contig(seqs: &mut Vec<(String, i64)>, index: &mut HashMap<String, usize>, c: &MafComp) {
    if !index.contains_key(&c.src) {
        index.insert(c.src.clone(), seqs.len());
        seqs.push((c.src.clone(), c.src_size as i64));
    }
}

/// Write a single two-sequence MAF block as a `.1aln` alignment record.
fn write_block(
    w: &mut pgr::libs::onepack::container::Writer,
    a: &MafComp,
    b: &MafComp,
    a_index: &HashMap<String, usize>,
    b_index: &HashMap<String, usize>,
    tspace: i64,
) -> Result<()> {
    // The CIGAR walks the `s` texts as given. `cigar_from_alignment(ref, qry)`
    // reports `qry` as the query, so pass `b` as ref and `a` as qry to make the
    // CIGAR a-vs-b (query = `a`, target = `b`), matching `paf to-1aln` and the
    // `.1aln` convention. When `b` is on the `-` strand its text is already the
    // reverse complement, matching the `.1aln` `R` marker.
    let ops = cigar_from_alignment(b.text.as_bytes(), a.text.as_bytes())?;
    if ops.is_empty() {
        bail!(
            "MAF block {}({}) vs {} has an empty alignment",
            a.src,
            a.strand,
            b.src
        );
    }

    let aread = *a_index
        .get(&a.src)
        .ok_or_else(|| anyhow!("MAF src {} not in skeleton", a.src))?;
    let abpos = a.start as i64;
    let aepos = (a.start + a.size) as i64;

    let bread = *b_index
        .get(&b.src)
        .ok_or_else(|| anyhow!("MAF src {} not in skeleton", b.src))?;
    let comp = b.strand == '-';
    let (bbpos, bepos) = if comp {
        let (s, e) = reverse_range_pair(b.start, b.start + b.size, b.src_size);
        (s as i64, e as i64)
    } else {
        (b.start as i64, (b.start + b.size) as i64)
    };

    // Sanity-check the CIGAR spans against the intervals.
    let a_span: i64 = ops.iter().map(|op| op.query_delta() as i64).sum();
    let b_span: i64 = ops.iter().map(|op| op.target_delta() as i64).sum();
    if a_span != aepos - abpos || b_span != bepos - bbpos {
        bail!(
            "MAF block {}({}) vs {}({}) CIGAR span {}:{} does not match intervals a[{}:{}] b[{}:{}]",
            a.src,
            a.strand,
            b.src,
            b.strand,
            a_span,
            b_span,
            abpos,
            aepos,
            bbpos,
            bepos
        );
    }

    // For a reverse-strand `b` side the stored path walks forward, so resample
    // the reversed CIGAR.
    let effective = if comp { reverse_cigar(&ops) } else { ops };
    let tp = cigar2tp(&effective, abpos, bbpos, tspace);

    write_aln_record(
        w,
        aread as i64,
        abpos,
        aepos,
        bread as i64,
        bbpos,
        bepos,
        comp,
        tp.diffs,
        &tp.tpoints,
        &tp.tdiffs,
    )
}
