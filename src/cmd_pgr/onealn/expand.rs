//! Shared read-side CLI layer for `pgr 1aln to-paf` / `to-psl`.
//!
//! Loads the source genomes referenced by a `.1aln` file and opens the file
//! (arg parsing + I/O setup). The record→PAF/PSL conversion itself lives in
//! `libs::onepack::expand`; this module only glues the clap args to it, per the
//! layering principle in `AGENTS.md`.

use anyhow::Context;
use std::io::Write;

use pgr::libs::onepack::expand::Genomes;
use pgr::libs::onepack::record::AlnFile;
use pgr::libs::paf::record::PafRecord;
use pgr::libs::pgi::build::read_fasta;

/// Open the `.1aln` file and load the two source genomes.
pub fn open_aln(args: &clap::ArgMatches) -> anyhow::Result<(AlnFile, Genomes)> {
    let infile = args
        .get_one::<String>("infile")
        .context("missing required argument: infile")?;
    let ref_seq = args
        .get_one::<String>("ref_seq")
        .context("missing required argument: --ref-seq")?;
    let qry_seq = args
        .get_one::<String>("query_seq")
        .context("missing required argument: --query-seq")?;
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    crate::cmd_pgr::args::ensure_outfile_distinct(
        outfile,
        [infile.as_str(), ref_seq.as_str(), qry_seq.as_str()],
    )?;

    let aln =
        AlnFile::open(infile).with_context(|| format!("Failed to open .1aln file {infile}"))?;
    let genomes = Genomes {
        ref_seqs: read_fasta(ref_seq)
            .with_context(|| format!("Failed to load reference FASTA {ref_seq}"))?,
        qry_seqs: read_fasta(qry_seq)
            .with_context(|| format!("Failed to load query FASTA {qry_seq}"))?,
    };
    Ok((aln, genomes))
}

/// Write a single PAF record to `out`.
pub fn write_paf<W: Write>(out: &mut W, rec: &PafRecord) -> anyhow::Result<()> {
    pgr::libs::paf::record::write_paf_record(out, rec)?;
    Ok(())
}
