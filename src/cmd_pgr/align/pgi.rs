//! `pgr align pgi` — pairwise genome alignment on the pgi k-mer pipeline.

use anyhow::Context;
use clap::parser::ValueSource;
use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};
use std::io::Read;
use std::path::{Path, PathBuf};

/// Build the clap subcommand for pgi.
pub fn make_subcommand() -> Command {
    Command::new("pgi")
        .about("Aligns two genomes or .pgi indexes into PSL blocks")
        .after_help(
            r###"
Merges the sorted k-mer streams of two genomes, chains the shared seeds in
anti-diagonal space, and emits one PSL block per chain. Block-level output is
meant to be chained by `pgr psl to_chain` / `pgr pl chainnet`.

Inputs may be genome sequences (FASTA, gzipped FASTA or .2bit) or .pgi indexes,
mixed freely:
* A genome sequence is indexed automatically: an index named like the input
  (e.g. ref.fa -> ref.pgi) is reused when present, otherwise one is built in a
  temporary directory (or next to the input with --keep-index). The sequence
  itself is then used to refine the chains.
* A .pgi index is used directly; --ref-seq/--query-seq may then supply the
  sequences for chain refinement, and are validated against the index.
With a single input the genome is aligned to itself (internal repeats and
haplotype-level homology, FastGA's self mode); exact self-identity hits are
dropped.

Notes:
* Both sides must use identical sampling parameters (k, syncmer, window).
* The query index is memory-mapped and must be a regular file ('stdin' and
  gzipped indexes are not supported).
* --k/--smer/--window apply only to genome-sequence inputs; .pgi inputs carry
  their parameters in the index header.
* K-mers occurring more than --freq times on either side are skipped.
* Chains shorter than --min-span on either axis are dropped.
* --ref-seq/--query-seq accept FASTA (.fa/.fa.gz) or .2bit files.

Examples:
1. Align two genomes directly (indexes built automatically):
   pgr align pgi ref.fa query.fa -o out.psl
2. Tune seed filtering and chaining:
   pgr align pgi ref.fa query.fa -f 20 -c 100 -s 2000 --band 64 -o out.psl
3. Reuse self-built indexes:
   pgr pgi build ref.fa -o ref.pgi
   pgr pgi build query.fa -o query.pgi
   pgr align pgi ref.pgi query.pgi --ref-seq ref.fa --query-seq query.fa -o out.psl
4. Stitch chains across small insertions:
   pgr align pgi ref.fa query.fa --merge-gap 10000 -o out.psl
5. Lower the partial-seed floor (tube defaults to FastGA's plen floor of 12):
   pgr align pgi ref.fa query.fa --min-shared 16 -o out.psl
6. Keep the automatically built indexes:
   pgr align pgi ref.fa query.fa --keep-index -o out.psl
"###,
        )
        .arg(
            Arg::new("ref")
                .index(1)
                .required(true)
                .help("Reference genome (FASTA/2bit) or .pgi index"),
        )
        .arg(
            Arg::new("query")
                .index(2)
                .help("Query genome (FASTA/2bit) or .pgi index; omit for self-alignment"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg_required())
        .arg(
            Arg::new("freq")
                .short('f')
                .long("freq")
                .default_value("10")
                .value_parser(value_parser!(u32))
                .help("Maximum k-mer frequency on either side to keep as seed"),
        )
        .arg(
            Arg::new("min_span")
                .short('c')
                .long("min-span")
                .default_value("85")
                .value_parser(value_parser!(u32))
                .help("Minimum per-axis seed span (bp) for a chain"),
        )
        .arg(
            Arg::new("max_gap")
                .short('s')
                .long("max-gap")
                .default_value("1000")
                .value_parser(value_parser!(u32))
                .help("Maximum bp gap between consecutive seeds in a chain"),
        )
        .arg(
            Arg::new("band")
                .long("band")
                .default_value("128")
                .value_parser(value_parser!(u32))
                .help("Diagonal band half-width (bp) around the chain mean"),
        )
        .arg(
            Arg::new("merge_gap")
                .long("merge-gap")
                .default_value("5000")
                .value_parser(value_parser!(u32))
                .help("Maximum gap (bp) between adjacent colinear chains to merge"),
        )
        .arg(
            Arg::new("min_shared")
                .long("min-shared")
                .value_parser(value_parser!(usize))
                .help("Minimum shared seed length (bp); default = k for greedy, 12 for tube"),
        )
        .arg(
            Arg::new("workflow")
                .long("workflow")
                .default_value("greedy")
                .value_parser(["greedy", "tube"])
                .help("Chaining workflow: greedy chains (default) or FastGA tubes"),
        )
        .arg(
            Arg::new("kmer")
                .short('k')
                .long("kmer")
                .default_value("40")
                .value_parser(value_parser!(usize))
                .help("k-mer size for automatic indexing (genome inputs only)"),
        )
        .arg(
            Arg::new("smer")
                .long("smer")
                .default_value("8")
                .value_parser(value_parser!(usize))
                .help("Syncmer s-mer length for automatic indexing"),
        )
        .arg(
            Arg::new("window")
                .long("window")
                .default_value("5")
                .value_parser(value_parser!(usize))
                .help("Syncmer window for automatic indexing"),
        )
        .arg(
            Arg::new("ref_seq").long("ref-seq").help(
                "Reference sequence file (FASTA or .2bit) for chain refinement of .pgi inputs",
            ),
        )
        .arg(
            Arg::new("query_seq")
                .long("query-seq")
                .help("Query sequence file (FASTA or .2bit) for chain refinement of .pgi inputs"),
        )
        .arg(
            Arg::new("keep_index")
                .long("keep-index")
                .action(ArgAction::SetTrue)
                .help("Keep automatically built indexes next to the genome inputs"),
        )
        .arg(crate::cmd_pgr::args::parallel_arg_with_default("8"))
}

/// One side of the alignment, resolved to an index path plus (for genome
/// inputs) the sequences themselves.
struct SideInput {
    index: String,
    seqs: Option<Vec<(String, Vec<u8>)>>,
}

/// Execute the align command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let ref_input = args.get_one::<String>("ref").unwrap();
    let query_input = args.get_one::<String>("query");
    let self_mode = query_input.is_none();
    let query_input = query_input.map(|s| s.as_str()).unwrap_or(ref_input);
    let outfile = args.get_one::<String>("outfile").unwrap();
    let params = pgr::libs::pgi::align::AlignParams {
        freq: *args.get_one::<u32>("freq").unwrap(),
        min_span: *args.get_one::<u32>("min_span").unwrap(),
        max_gap: *args.get_one::<u32>("max_gap").unwrap(),
        band: *args.get_one::<u32>("band").unwrap(),
        merge_gap: *args.get_one::<u32>("merge_gap").unwrap(),
        min_shared: args.get_one::<usize>("min_shared").copied(),
        workflow: match args.get_one::<String>("workflow").unwrap().as_str() {
            "tube" => pgr::libs::pgi::align::Workflow::Tube,
            _ => pgr::libs::pgi::align::Workflow::Greedy,
        },
    };
    let keep = args.get_flag("keep_index");
    let mut tmp: Option<tempfile::TempDir> = None;

    let SideInput {
        index: ref_index,
        seqs: ref_side_seqs,
    } = resolve_side(args, ref_input, "reference", &mut tmp, keep)?;
    let SideInput {
        index: query_index,
        seqs: query_side_seqs,
    } = if self_mode {
        // Self-alignment resolves the same input once and reuses it on both
        // sides (the sequence copy is bounded by the input size).
        SideInput {
            index: ref_index.clone(),
            seqs: ref_side_seqs.clone(),
        }
    } else {
        resolve_side(args, query_input, "query", &mut tmp, keep)?
    };
    let _tmp_guard = tmp;

    // The reference index is consumed as a stream by the merge; the query
    // index is memory-mapped (FastGA's GIX model) and decoded on demand, so
    // neither index is materialized in full.
    let mut r1 = pgr::reader(&ref_index)?;
    let mut a = pgr::libs::pgi::PgiStream::open(&mut r1)?;
    let b = pgr::libs::pgi::PgiMmap::open(std::path::Path::new(&query_index))?;

    // Extension sequences come from genome inputs directly (validated against
    // the index, which matters when a sibling index was reused) or from
    // --ref-seq/--query-seq for .pgi inputs (validated the same way).
    let ref_seqs = resolve_seqs(
        args,
        ref_side_seqs,
        a.header().contigs.as_slice(),
        "reference",
        "ref_seq",
    )?;
    let query_seqs = resolve_seqs(args, query_side_seqs, b.contigs(), "query", "query_seq")?;
    if ref_seqs.is_empty() != query_seqs.is_empty() {
        anyhow::bail!(
            "extension sequences are needed for both sides (genome inputs, or \
             --ref-seq/--query-seq for .pgi inputs)"
        );
    }

    let parallel = *args.get_one::<usize>("parallel").unwrap();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(parallel)
        .build()
        .context("building align thread pool")?;
    let psls = pool.install(|| -> anyhow::Result<Vec<pgr::libs::fmt::psl::Psl>> {
        if ref_seqs.is_empty() {
            pgr::libs::pgi::align::align_to_psl_streaming(&mut a, &b, &params, self_mode)
        } else {
            pgr::libs::pgi::align::align_to_psl_ext_streaming(
                a,
                b,
                &params,
                &ref_seqs,
                &query_seqs,
                self_mode,
            )
        }
    })?;
    let mut writer = pgr::writer(outfile)?;
    for p in &psls {
        p.write_to(&mut writer)?;
    }
    log::info!(
        "wrote {} PSL blocks (freq={}, min-span={}, max-gap={}, band={}) to {}",
        psls.len(),
        params.freq,
        params.min_span,
        params.max_gap,
        params.band,
        outfile
    );
    Ok(())
}

/// Resolve one input to an index path: .pgi inputs are used directly; genome
/// inputs reuse a same-named sibling .pgi or build one (temporarily, or next
/// to the input with `--keep-index`).
fn resolve_side(
    args: &ArgMatches,
    input: &str,
    label: &str,
    tmp: &mut Option<tempfile::TempDir>,
    keep: bool,
) -> anyhow::Result<SideInput> {
    if is_pgi_input(input) {
        return Ok(SideInput {
            index: input.to_string(),
            seqs: None,
        });
    }

    let input_path = Path::new(input);
    let cached = sibling_pgi_path(input_path);
    if cached.exists() {
        let (ck, cs, cw) = read_index_params(&cached)?;
        let explicit = |name: &str| args.value_source(name) == Some(ValueSource::CommandLine);
        let k = *args.get_one::<usize>("kmer").unwrap();
        let smer = *args.get_one::<usize>("smer").unwrap();
        let window = *args.get_one::<usize>("window").unwrap();
        if explicit("kmer") && k != ck {
            anyhow::bail!(
                "--kmer {k} conflicts with the cached index {} (k={ck})",
                cached.display()
            );
        }
        if explicit("smer") && smer != cs {
            anyhow::bail!(
                "--smer {smer} conflicts with the cached index {} (smer={cs})",
                cached.display()
            );
        }
        if explicit("window") && window != cw {
            anyhow::bail!(
                "--window {window} conflicts with the cached index {} (window={cw})",
                cached.display()
            );
        }
        log::info!("reusing {label} index {}", cached.display());
        let seqs = read_seqs(input)?;
        return Ok(SideInput {
            index: cached.display().to_string(),
            seqs: Some(seqs),
        });
    }

    let k = *args.get_one::<usize>("kmer").unwrap();
    let smer = *args.get_one::<usize>("smer").unwrap();
    let window = *args.get_one::<usize>("window").unwrap();
    let seqs = read_seqs(input)?;
    let idx = pgr::libs::pgi::build::build_from_seqs(seqs.clone(), k, smer, window, false)?;
    let out = if keep {
        cached.clone()
    } else {
        let dir = tmp.get_or_insert_with(|| {
            tempfile::TempDir::new().expect("creating temporary index directory")
        });
        dir.path().join(format!("{label}.pgi"))
    };
    let mut w = std::fs::File::create(&out)?;
    idx.write(&mut w)?;
    log::info!(
        "built {label} index {} (k={k}, syncmer {smer}/{window})",
        out.display()
    );
    Ok(SideInput {
        index: out.display().to_string(),
        seqs: Some(seqs),
    })
}

/// Resolve the extension sequences for one side, validating them against the
/// index contig table.
fn resolve_seqs(
    args: &ArgMatches,
    side_seqs: Option<Vec<(String, Vec<u8>)>>,
    header: &[(String, u64)],
    label: &str,
    arg_name: &str,
) -> anyhow::Result<Vec<(String, Vec<u8>)>> {
    match (side_seqs, args.get_one::<String>(arg_name)) {
        (Some(_), Some(_)) => anyhow::bail!(
            "{label} input is a genome sequence; --{} applies only to .pgi inputs",
            arg_name.replace('_', "-")
        ),
        (Some(seqs), None) => {
            validate_contigs(header, &seqs, label)?;
            Ok(seqs)
        }
        (None, Some(path)) => {
            let seqs = read_seqs(path)?;
            validate_contigs(header, &seqs, label)?;
            Ok(seqs)
        }
        (None, None) => Ok(Vec::new()),
    }
}

/// A .pgi input is detected by its magic; a missing file falls back to the
/// `.pgi` extension so the error stays on the index path.
fn is_pgi_input(path: &str) -> bool {
    let p = Path::new(path);
    if let Ok(mut f) = std::fs::File::open(p) {
        let mut magic = [0u8; 4];
        if f.read_exact(&mut magic).is_ok() {
            return &magic == pgr::libs::pgi::PGI_MAGIC;
        }
    }
    p.extension().and_then(|e| e.to_str()) == Some("pgi")
}

/// Sibling index path for a genome input: ref.fa / ref.fa.gz / ref.2bit all
/// map to ref.pgi.
fn sibling_pgi_path(path: &Path) -> PathBuf {
    let mut p = path.to_path_buf();
    if p.extension().and_then(|e| e.to_str()) == Some("gz") {
        p.set_extension("");
    }
    p.set_extension("pgi");
    p
}

/// Sampling parameters stored in a .pgi header.
fn read_index_params(path: &Path) -> anyhow::Result<(usize, usize, usize)> {
    let mut r = pgr::reader(path.to_str().context("index path utf8")?)?;
    let stream = pgr::libs::pgi::PgiStream::open(&mut r)?;
    let h = stream.header();
    Ok((h.k, h.smer, h.window))
}

/// The contig table of a .pgi index must match the extension sequences.
fn validate_contigs(
    header: &[(String, u64)],
    seqs: &[(String, Vec<u8>)],
    label: &str,
) -> anyhow::Result<()> {
    anyhow::ensure!(
        header.len() == seqs.len(),
        "{label} sequence count {} does not match index contig count {}",
        seqs.len(),
        header.len()
    );
    for ((name, len), (sname, s)) in header.iter().zip(seqs) {
        anyhow::ensure!(
            name == sname,
            "{label} contig name mismatch: index '{name}' vs sequence '{sname}'"
        );
        anyhow::ensure!(
            *len == s.len() as u64,
            "{label} contig '{name}' length mismatch: index {len} vs sequence {}",
            s.len()
        );
    }
    Ok(())
}

/// Read all sequences from a FASTA (plain or gzipped) or .2bit file.
fn read_seqs(path: &str) -> anyhow::Result<Vec<(String, Vec<u8>)>> {
    let is_2bit = std::path::Path::new(path)
        .extension()
        .and_then(|e| e.to_str())
        == Some("2bit");
    if is_2bit {
        pgr::libs::pgi::build::read_2bit(path, false)
    } else {
        pgr::libs::pgi::build::read_fasta(path)
    }
}
