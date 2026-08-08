//! Create a new pbit archive from a reference FASTA and sample FASTA files.

use anyhow::Context;
use clap::{ArgMatches, Command};
use pgr::libs::pbit::compressor::Compressor;
use pgr::libs::pbit::format::MAX_PACKED_SIZE;

/// Build the clap subcommand for create.
pub fn make_subcommand() -> Command {
    Command::new("create")
        .about("Creates a pbit archive from a reference FASTA and sample FASTA files")
        .after_help(
            r###"
This command creates a new pbit archive. The reference FASTA is stored as
standard 2bit records; each sample FASTA is encoded against the matching
reference segment and stored as delta entries.

PAF is mandatory (`--paf`, or the TSV 3rd column): segments covered by PAF
alignments are CIGAR-encoded (pure-match segments become zero-cost Identity
references); uncovered segments fall back to LZ-diff/Raw. An empty PAF file
disables CIGAR encoding (all segments use LZ-diff/Raw).

Notes:
* Sample names are derived from the input FASTA basenames (use `--name` to
  override with a TSV file of `name<TAB>path[<TAB>paf_path][<TAB>ref_name]` lines)
* Reference and sample FASTA files may be plain text or gzipped (.gz)
* Sequences with no matching reference content are stored verbatim (Raw
  deltas), so the archive is lossless for ACGTN input; degenerate IUPAC
  codes are the only accepted loss (see below)
* Only ACGTN characters are supported; IUPAC degenerate codes (R, Y, S, W,
  K, M, B, D, H, V) are lossily mapped to N
* `--paf` files are paired with `-i` files by order; `--name` and `--paf`
  are mutually exclusive (use the TSV's 3rd column for PAF, which is required)

Examples:
1. Create a pbit archive with one sample:
   pgr pbit create -r ref.fa -i sample1.fa -p sample1.paf -o out.pbit

2. Create with multiple samples:
   pgr pbit create -r ref.fa -i s1.fa -i s2.fa -i s3.fa -o out.pbit

3. Custom segment size and k-mer length:
   pgr pbit create -r ref.fa -i sample.fa -o out.pbit -s 8192 -k 15

4. Provide sample names via a TSV file:
   pgr pbit create -r ref.fa --name samples.tsv -o out.pbit

5. CIGAR-driven encoding with PAF:
   pgr pbit create -r ref.fa -i sample.fa -p sample.paf -o out.pbit

6. Multiple references (samples route to reference 0 by default, or via the
   TSV's 4th column):
   pgr pbit create -r ref1.fa -r ref2.fa --name samples.tsv -o out.pbit
"###,
        )
        .arg(crate::cmd_pgr::args::pbit_ref_arg())
        .arg(crate::cmd_pgr::args::pbit_infiles_arg())
        .arg(crate::cmd_pgr::args::outfile_arg_required())
        .arg(crate::cmd_pgr::args::pbit_segment_size_arg())
        .arg(crate::cmd_pgr::args::pbit_kmer_len_arg())
        .arg(crate::cmd_pgr::args::pbit_min_match_len_arg())
        .arg(crate::cmd_pgr::args::pbit_name_arg())
        .arg(crate::cmd_pgr::args::pbit_paf_arg())
}

/// Execute the create command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let ref_fastas: Vec<&str> = args
        .get_many::<String>("ref")
        .context("missing required argument: --ref")?
        .map(|s| s.as_str())
        .collect();
    let outfile = args
        .get_one::<String>("outfile")
        .context("missing required argument: --outfile")?;
    let segment_size = *args
        .get_one::<usize>("segment_size")
        .context("missing --segment-size")?;
    let kmer_len = *args
        .get_one::<usize>("kmer_len")
        .context("missing --kmer-len")?;
    let min_match_len = *args
        .get_one::<u32>("min_match_len")
        .context("missing --min-match-len")?;

    anyhow::ensure!(segment_size > 0, "segment-size must be positive");
    anyhow::ensure!(
        segment_size <= i32::MAX as usize,
        "segment-size must not exceed {}",
        i32::MAX
    );
    anyhow::ensure!(kmer_len > 0, "kmer-len must be positive");
    anyhow::ensure!(min_match_len > 0, "min-match-len must be positive");
    // A match/k-mer cannot span more than one reference segment, so lengths
    // greater than segment-size are meaningless. More importantly, min_match_len
    // drives `LzDiff::prepare`'s `reference.resize(len + key_len)` padding: an
    // unbounded value (e.g. `-l 4294967295`) would trigger a multi-GB allocation
    // per segment. Bound both to keep the CLI sane.
    anyhow::ensure!(
        kmer_len <= segment_size,
        "kmer-len ({}) must not exceed segment-size ({})",
        kmer_len,
        segment_size
    );
    anyhow::ensure!(
        min_match_len as usize <= segment_size,
        "min-match-len ({}) must not exceed segment-size ({})",
        min_match_len,
        segment_size
    );
    // `min_match_len` drives `LzDiff::prepare`'s `reference.resize(len + key_len)`
    // padding (`key_len ≈ min_match_len`), so an unbounded value would force a
    // multi-GB allocation per decoded segment. Apply the same absolute cap the
    // decompressor enforces (`MAX_PACKED_SIZE`) so `create` never produces an
    // archive that `stat` / `range` / `to-fa` would later reject as invalid.
    anyhow::ensure!(
        min_match_len as usize <= MAX_PACKED_SIZE,
        "min-match-len ({}) must not exceed the per-segment bound ({})",
        min_match_len,
        MAX_PACKED_SIZE
    );

    // Guard against -o truncating an input file before it is read (e.g. a
    // reference or sample FASTA, PAF, or the --name TSV). create_multi opens
    // the output with File::create, which truncates, BEFORE reading inputs.
    let samples = super::collect_samples_from_args(args)?;
    let mut inputs: Vec<&str> = ref_fastas.clone();
    for (_, path, paf_opt, _) in &samples {
        inputs.push(path.as_str());
        if let Some(paf) = paf_opt {
            inputs.push(paf.as_str());
        }
    }
    if let Some(name_tsv) = args.get_one::<String>("name") {
        inputs.push(name_tsv.as_str());
    }
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, inputs)?;

    let mut comp =
        Compressor::create_multi(outfile, &ref_fastas, segment_size, kmer_len, min_match_len)
            .with_context(|| format!("failed to create pbit archive: {}", outfile))?;

    let mut cmd_line = format!(
        "pgr pbit create -r {} -o {} -s {} -k {} -l {}",
        ref_fastas.join(" -r "),
        outfile,
        segment_size,
        kmer_len,
        min_match_len
    );
    // Record the sample inputs for provenance (the gathered sample names and
    // paths are not re-parsed — they document which FASTA/PAF fed the archive).
    for (name, path, paf_opt, ref_spec) in &samples {
        cmd_line.push_str(&format!(" -i {}:{}", name, path));
        if let Some(paf) = paf_opt {
            cmd_line.push_str(&format!(" -p {}", paf));
        }
        if let Some(ref_spec) = ref_spec {
            cmd_line.push_str(&format!(" @ref {}", ref_spec));
        }
    }
    comp.set_cmd_line(&cmd_line);

    for (name, path, paf_opt, ref_spec) in &samples {
        let ref_id = super::resolve_ref_id(ref_spec.as_deref(), &ref_fastas)?;
        comp.set_cur_ref_id(ref_id);
        match paf_opt {
            Some(paf) => comp
                .append_sample_with_paf(name, path, paf)
                .with_context(|| format!("failed to append sample '{}' with PAF", name))?,
            None => comp
                .append_sample(name, path)
                .with_context(|| format!("failed to append sample '{}'", name))?,
        }
    }
    comp.finish().context("failed to finalize pbit archive")?;

    Ok(())
}
