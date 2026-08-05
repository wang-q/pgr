//! `pgr paf to-1aln`: convert a PAF file with `cg:Z` CIGAR into FastGA `.1aln`.
//!
//! This is the write-side mirror of `pgr 1aln to-paf`. It resamples each
//! record's X-CIGAR into trace points (`cigar2tp`) and writes the GDB skeleton
//! plus alignment records into a ONEcode `.1aln` container, following the
//! design doc `notes/design/1aln.md` §6.1.1 / §7.9 (FastGA `PAFtoALN`).
//!
//! Axis convention (matching `PAFtoALN` and `pgr 1aln to-paf`): the PAF query
//! is the `.1aln` `a` side and the PAF target is the `b` side. A `-` strand
//! marks the `b` side reverse-complemented (`R` line), and the target
//! coordinates are reversed onto the forward source.

use anyhow::{anyhow, bail, Context, Result};
use clap::{ArgMatches, Command};
use std::collections::HashMap;
use std::io::BufReader;

use pgr::libs::onepack::write::{
    cigar2tp, open_aln_writer, write_aln_record, write_skeleton_contigs, write_tspace,
};
use pgr::libs::paf::cigar::{extract_cigar, reverse_cigar};
use pgr::libs::paf::parser::parse_paf;
use pgr::libs::paf::record::PafRecord;

/// Build the clap subcommand for to-1aln.
pub fn make_subcommand() -> Command {
    Command::new("to-1aln")
        .about("Converts a PAF file with cg:Z CIGAR to FastGA .1aln format")
        .after_help(
            r###"
Resamples each PAF record's `cg:Z` X-CIGAR into trace points (`cigar2tp`) and
writes them, along with the GDB skeleton (one contig per PAF sequence name),
into a FastGA `.1aln` (ONEcode trace-point) file.

The PAF query becomes the `.1aln` `a` side and the PAF target the `b` side. A
`-` strand marks the `b` side reverse-complemented (`R` line); the target
coordinates are reversed onto the forward source and the CIGAR is reversed.
This is the write-side mirror of `pgr 1aln to-paf`.

Notes:
* Every PAF record must carry a `cg:Z` tag with `=`/`X`/`I`/`D` ops (M is
  treated as `=`). Records without a CIGAR are an error.
* The ONEcode container is binary and requires a real output file (not stdout).
* The skeleton is built from the PAF sequence names and lengths; no source
  FASTA is required.

Examples:
1. Convert a PAF with CIGAR to .1aln:
   pgr paf to-1aln align.paf -o align.1aln
2. Use a custom trace-point spacing:
   pgr paf to-1aln align.paf -o align.1aln --tspace 50

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input PAF file",
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

    let file =
        std::fs::File::open(infile).with_context(|| format!("Failed to open PAF file {infile}"))?;
    let records = parse_paf(BufReader::new(file))
        .with_context(|| format!("Failed to parse PAF file {infile}"))?;
    if records.is_empty() {
        bail!("no PAF records found in {infile}");
    }

    let mut writer = open_aln_writer(outfile)?;
    // Reference the source PAF file (db1, count 1) before the header flushes.
    writer.add_reference(infile, 1);
    write_tspace(&mut writer, tspace)?;

    // Build the GDB skeleton from the PAF sequence names/lengths. The query
    // side is the `a` contigs, the target side the `b` contigs.
    let (a_seqs, a_index) =
        collect_sequences(&records, |r| (r.query_name.clone(), r.query_length as i64));
    let (b_seqs, b_index) = collect_sequences(&records, |r| {
        (r.target_name.clone(), r.target_length as i64)
    });

    write_skeleton_contigs(&mut writer, &a_seqs)?;
    write_skeleton_contigs(&mut writer, &b_seqs)?;

    for rec in &records {
        write_record(&mut writer, rec, &a_index, &b_index, tspace)?;
    }

    writer.close()?;
    Ok(())
}

/// Collect the distinct sequence names in first-seen order.
///
/// Returns the `(name, length)` list (order = contig index) and a name → index
/// map. `name_fn` picks the query or target `(name, length)` from a record.
fn collect_sequences<F>(
    records: &[PafRecord],
    name_fn: F,
) -> (Vec<(String, i64)>, HashMap<String, usize>)
where
    F: Fn(&PafRecord) -> (String, i64),
{
    let mut index = HashMap::new();
    let mut seqs: Vec<(String, i64)> = Vec::new();
    for rec in records {
        let (name, len) = name_fn(rec);
        if !index.contains_key(&name) {
            index.insert(name.clone(), seqs.len());
            seqs.push((name, len));
        }
    }
    (seqs, index)
}

/// Write a single PAF record as a `.1aln` alignment record.
///
/// Query is the `a` side (always forward); target is the `b` side. On `-`
/// strand the target coordinates are reversed onto the forward source and the
/// CIGAR is reversed so the stored path walks forward on both axes.
fn write_record(
    w: &mut pgr::libs::onepack::container::Writer,
    rec: &PafRecord,
    a_index: &HashMap<String, usize>,
    b_index: &HashMap<String, usize>,
    tspace: i64,
) -> Result<()> {
    let ops = extract_cigar(&rec.tags)?;
    if ops.is_empty() {
        bail!(
            "PAF record {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{} is missing a cg:Z CIGAR",
            rec.query_name,
            rec.query_length,
            rec.query_start,
            rec.query_end,
            rec.strand,
            rec.target_name,
            rec.target_length,
            rec.target_start,
            rec.target_end,
            rec.matches,
            rec.block_length,
            rec.mapq
        );
    }

    let aread = *a_index
        .get(&rec.query_name)
        .ok_or_else(|| anyhow!("query {} not in skeleton", rec.query_name))?;
    let abpos = rec.query_start as i64;
    let aepos = rec.query_end as i64;

    let bread = *b_index
        .get(&rec.target_name)
        .ok_or_else(|| anyhow!("target {} not in skeleton", rec.target_name))?;
    let comp = rec.strand == '-';
    let (bbpos, bepos) = if comp {
        // Reverse the target interval onto the forward source.
        let bbpos = rec.target_length as i64 - rec.target_end as i64;
        let bepos = rec.target_length as i64 - rec.target_start as i64;
        (bbpos, bepos)
    } else {
        (rec.target_start as i64, rec.target_end as i64)
    };

    // Sanity-check the CIGAR spans against the intervals.
    let a_span: i64 = ops.iter().map(|op| op.query_delta() as i64).sum();
    let b_span: i64 = ops.iter().map(|op| op.target_delta() as i64).sum();
    if a_span != aepos - abpos || b_span != bepos - bbpos {
        bail!(
            "PAF record {}\t{}\tcig:Z:{}\t span {}:{} does not match intervals a[{}:{}] b[{}:{}]",
            rec.query_name,
            rec.target_name,
            pgr::libs::paf::cigar::format_cigar(&ops),
            a_span,
            b_span,
            abpos,
            aepos,
            bbpos,
            bepos
        );
    }

    // For a reverse-strand record the stored path walks the target forward, so
    // resample the reversed CIGAR.
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
