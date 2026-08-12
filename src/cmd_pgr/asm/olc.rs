use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use pgr::libs::asm::assemble::{assemble_unitigs_buf, AssembleOptions};
use pgr::libs::olc::consensus::consensus;
use pgr::libs::olc::layout::build_layouts;
use pgr::libs::olc::overlap::{find_overlaps, OverlapOptions, Unitig};
use std::io::Write;
use std::path::Path;

/// Build the clap subcommand for olc.
pub fn make_subcommand() -> Command {
    Command::new("olc")
        .about("Assembles reads into contigs via multi-k unitig OLC")
        .after_help(
            r###"
Runs the full OLC pipeline in memory: for every k in --kmer the reads are
assembled into maximal unitigs (`pgr asm unitig` semantics), all unitigs are
pooled as pseudo-reads, exact overlaps are found (`pgr asm ovlp`), layouts
are built greedily (`pgr asm layout`), and each layout is stitched into a
consensus contig (`pgr asm cns`). See notes/design/olc.md.

Unitigs are named `k<k>:unitig_<id>` so the per-k sets stay distinguishable
and reproducible. Overlaps are exact (error-free unitigs), layouts stop at
ambiguous junctions and non-reciprocal edges, and no bubble heuristics are
applied.

Notes:
* Input is 1 interleaved file or 2 paired files (same as `pgr asm unitig`)
* --keep-dir writes the intermediate unitigs/overlap/layout files for
  debugging or for re-running the stage commands separately
* Output contigs are written longest-first with `>contig_<id>,len=...,cov=...`
  headers, 70-column wrapped

Examples:
1. Assemble a small metagenome with three k values:
   pgr asm olc reads.fq.gz -o contigs.fa --kmer 21,51,81
2. Keep the intermediates and raise the minimum contig length:
   pgr asm olc R1.fq.gz R2.fq.gz -o contigs.fa \
       --kmer 21,51,81 --min-contig-len 1000 --keep-dir stage/
"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg_with_numargs(
            "Input reads: 1 interleaved or 2 paired (R1, R2)",
            1..=2,
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("kmer")
                .long("kmer")
                .short('k')
                .num_args(1)
                .default_value("21,51,81")
                .help("Comma-separated k-mer lengths for the unitig sets"),
        )
        .arg(
            Arg::new("min_count_seed")
                .long("min-count-seed")
                .num_args(1)
                .default_value("3")
                .value_parser(value_parser!(usize))
                .help("Solid k-mer count threshold for unitig assembly"),
        )
        .arg(
            Arg::new("overlap_k")
                .long("overlap-k")
                .num_args(1)
                .default_value("17")
                .value_parser(value_parser!(usize))
                .help("Seed k-mer length for overlap detection"),
        )
        .arg(
            Arg::new("min_overlap")
                .long("min-overlap")
                .num_args(1)
                .default_value("34")
                .value_parser(value_parser!(usize))
                .help("Minimum accepted overlap length in bases"),
        )
        .arg(
            Arg::new("min_contig_len")
                .long("min-contig-len")
                .num_args(1)
                .default_value("500")
                .value_parser(value_parser!(usize))
                .help("Minimum output contig length in bases"),
        )
        .arg(
            Arg::new("keep_dir")
                .long("keep-dir")
                .num_args(1)
                .help("Directory for intermediate unitigs/ovlp/layout files"),
        )
}

/// Execute the olc command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let ks: Vec<usize> = args
        .get_one::<String>("kmer")
        .unwrap()
        .split(',')
        .map(|s| {
            s.trim()
                .parse::<usize>()
                .with_context(|| format!("invalid --kmer value: {s}"))
        })
        .collect::<anyhow::Result<_>>()?;
    anyhow::ensure!(!ks.is_empty(), "at least one k-mer length is required");
    let min_count_seed = *args.get_one::<usize>("min_count_seed").unwrap();
    let seed_k = *args.get_one::<usize>("overlap_k").unwrap();
    let min_overlap = *args.get_one::<usize>("min_overlap").unwrap();
    let min_contig_len = *args.get_one::<usize>("min_contig_len").unwrap();
    let keep_dir = args.get_one::<String>("keep_dir");
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    // S0: unitigs per k.
    let mut unitigs = Vec::new();
    for &k in &ks {
        let opts = AssembleOptions {
            k,
            min_count_seed,
            ..AssembleOptions::default()
        };
        let bufs = assemble_unitigs_buf(&infiles, &opts)?;
        for (id, bases) in bufs {
            unitigs.push(Unitig {
                name: format!("k{k}:unitig_{id}"),
                seq: bases,
            });
        }
    }
    if let Some(dir) = keep_dir {
        std::fs::create_dir_all(dir)
            .with_context(|| format!("failed to create --keep-dir {dir}"))?;
        dump_unitigs(dir, &unitigs)?;
    }

    // S1: exact overlaps.
    let overlaps = find_overlaps(
        &unitigs,
        &OverlapOptions {
            seed_k,
            min_overlap,
        },
    )?;
    if let Some(dir) = keep_dir {
        dump_paf(dir, &unitigs, &overlaps)?;
    }

    // S2: greedy layouts.
    let layouts = build_layouts(&unitigs, &overlaps);
    if let Some(dir) = keep_dir {
        dump_layouts(dir, &unitigs, &layouts)?;
    }

    // S3: consensus stitch.
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

/// Writes the pooled unitigs as one FASTA (k-tagged names).
fn dump_unitigs(dir: &str, unitigs: &[Unitig]) -> anyhow::Result<()> {
    let path = Path::new(dir).join("unitigs.fa");
    let path = path.to_str().unwrap();
    let mut out = pgr::libs::io::writer(path).with_context(|| format!("failed to open {path}"))?;
    for u in unitigs {
        writeln!(out, ">{}", u.name)?;
        for chunk in u.seq.chunks(70) {
            out.write_all(chunk)?;
            out.write_all(b"\n")?;
        }
    }
    out.flush()?;
    Ok(())
}

/// Writes the overlap PAF.
fn dump_paf(
    dir: &str,
    unitigs: &[Unitig],
    overlaps: &[pgr::libs::olc::overlap::Overlap],
) -> anyhow::Result<()> {
    let path = Path::new(dir).join("ovlp.paf");
    let path = path.to_str().unwrap();
    let mut out = pgr::libs::io::writer(path).with_context(|| format!("failed to open {path}"))?;
    for ov in overlaps {
        let rec = super::common::to_paf(ov, unitigs);
        pgr::libs::paf::record::write_paf_record(&mut out, &rec)?;
    }
    out.flush()?;
    Ok(())
}

/// Writes the layout TSV.
fn dump_layouts(
    dir: &str,
    unitigs: &[Unitig],
    layouts: &[pgr::libs::olc::layout::Layout],
) -> anyhow::Result<()> {
    let path = Path::new(dir).join("layout.tsv");
    let path = path.to_str().unwrap();
    let mut out = pgr::libs::io::writer(path).with_context(|| format!("failed to open {path}"))?;
    super::common::write_layout_tsv(&mut out, unitigs, layouts)?;
    out.flush()?;
    Ok(())
}
