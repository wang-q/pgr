use anyhow::Context;
use clap::{ArgMatches, Command};
use pgr::libs::onepack::record::AlnFile;

/// Build the clap subcommand for stat.
pub fn make_subcommand() -> Command {
    Command::new("stat")
        .about("Reports header and per-record statistics from a .1aln file")
        .after_help(
            r###"
Reads a FastGA `.1aln` (ONEcode trace-point) file and prints a TSV report
(key<TAB>value) of the header and record statistics. No source sequences are
needed.

Reported fields:
* tspace  - trace point spacing (the `t` line)
* records - number of `A` alignment records (from the footer counts)
* skeletons - number of GDB skeleton objects
* scaffolds / contigs - across all skeletons
* refs - number of reference (`<`) entries, each with count 1/2/3
* provenance - program/version/command that produced the file

Notes:
* Reads a single .1aln file; does not support gzip or stdin (the ONEcode
  container requires random access to the footer offset at EOF).

Examples:
1. Report header stats:
   pgr 1aln stat mg1655-sakai.1aln -o stats.tsv

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input .1aln file",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the stat command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args
        .get_one::<String>("infile")
        .context("missing required argument: infile")?;
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, [infile.as_str()])?;

    let aln =
        AlnFile::open(infile).with_context(|| format!("Failed to open .1aln file {infile}"))?;
    let mut records = 0i64;
    let mut max_points = 0usize;
    let mut total_diffs = 0i64;
    let mut aln = aln;
    while let Some(rec) = aln.next_record()? {
        records += 1;
        max_points = max_points.max(rec.num_points());
        total_diffs += rec.diffs;
    }
    let n_scaffolds: usize = aln.skeletons.iter().map(|s| s.scaffolds.len()).sum();
    let n_contigs: usize = aln.skeletons.iter().map(|s| s.contigs.len()).sum();

    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {outfile}"))?;
    writeln_report(&mut writer, "tspace", &aln.tspace)?;
    writeln_report(&mut writer, "records", &records)?;
    writeln_report(&mut writer, "max_trace_points", &max_points)?;
    writeln_report(&mut writer, "total_diffs", &total_diffs)?;
    writeln_report(&mut writer, "skeletons", &aln.skeletons.len())?;
    writeln_report(&mut writer, "scaffolds", &n_scaffolds)?;
    writeln_report(&mut writer, "contigs", &n_contigs)?;
    writeln_report(&mut writer, "refs", &aln.references().len())?;
    for (i, r) in aln.references().iter().enumerate() {
        writeln_report(&mut writer, &format!("ref.{i}.filename"), &r.filename)?;
        writeln_report(&mut writer, &format!("ref.{i}.count"), &r.count)?;
    }
    for (i, p) in aln.provenance().iter().enumerate() {
        writeln_report(&mut writer, &format!("prov.{i}.program"), &p.program)?;
        writeln_report(&mut writer, &format!("prov.{i}.version"), &p.version)?;
    }
    Ok(())
}

fn writeln_report<W, T>(w: &mut W, key: &str, value: &T) -> anyhow::Result<()>
where
    W: std::io::Write,
    T: std::fmt::Display,
{
    writeln!(w, "{key}\t{value}")?;
    Ok(())
}
