//! `pgr align pgi` — pairwise genome alignment on the pgi k-mer pipeline.

use anyhow::Context;
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
meant to be chained by `pgr psl to-chain` / `pgr pl chainnet`.

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
dropped. `--self` states the same explicitly and accepts the reference input
again as the query.

Notes:
* Both sides must use identical sampling parameters (k, syncmer, window).
* The query index is memory-mapped and must be a regular file ('stdin' and
  gzipped indexes are not supported).
* --kmer/--smer/--window apply only to genome-sequence inputs; .pgi inputs carry
  their parameters in the index header.
* K-mers occurring at least --freq times on either side are skipped.
* Without extension sequences (a .pgi pair without --ref-seq/--query-seq)
  each tube is emitted as one geometric block; with sequences the tubes are
  refined by FastGA's mid-line wave into scored multi-block records.
* --ref-seq/--query-seq accept FASTA (.fa/.fa.gz) or .2bit files.

Examples:
1. Align two genomes directly (indexes built automatically):
   pgr align pgi ref.fa query.fa -o out.psl
2. Tune seed filtering and chaining:
   pgr align pgi ref.fa query.fa -f 20 -o out.psl
3. Reuse self-built indexes:
   pgr pgi build ref.fa -o ref.pgi
   pgr pgi build query.fa -o query.pgi
   pgr align pgi ref.pgi query.pgi --ref-seq ref.fa --query-seq query.fa -o out.psl
4. Chain without extension sequences (geometric blocks):
   pgr align pgi ref.pgi query.pgi -o out.psl
5. Lower the partial-seed floor (default is FastGA's plen floor of 12):
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
        .arg(
            Arg::new("self_align")
                .long("self")
                .action(clap::ArgAction::SetTrue)
                .help("Self-alignment (query omitted or the same input as the reference)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg_required())
        .arg(
            Arg::new("freq")
                .short('f')
                .long("freq")
                .default_value("10")
                .value_parser(value_parser!(u32))
                .help("K-mers occurring at least this often on either side are skipped as seeds"),
        )
        .arg(
            Arg::new("min_shared")
                .long("min-shared")
                .value_parser(value_parser!(usize))
                .help("Minimum shared seed length (bp); default = FastGA's plen floor (12)"),
        )
        .arg(
            Arg::new("kmer")
                .short('k')
                .long("kmer")
                .default_value("40")
                .value_parser(value_parser!(usize))
                .help("k-mer size for automatic indexing (genome inputs only, must be <= 64)"),
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
    let is_self = args.get_flag("self_align");
    if is_self {
        if let Some(q) = query_input {
            anyhow::ensure!(
                q == ref_input,
                "--self expects the query to be the same input as the reference \
                 (omit the query or pass the same file)"
            );
        }
    }
    let self_mode = is_self || query_input.is_none();
    let query_input = query_input.map(|s| s.as_str()).unwrap_or(ref_input);
    let outfile = args.get_one::<String>("outfile").unwrap();
    let mut inputs: Vec<String> = vec![ref_input.to_string(), query_input.to_string()];
    // Also protect the sibling index each genome input maps to (`ref.fa` ->
    // `ref.pgi`, `ref.fa.gz` -> `ref.fa.pgi`): `-o ref.pgi` must not silently
    // overwrite the index the command creates or reuses, which would corrupt
    // it and break the next run with a confusing "reading header" error.
    for s in [ref_input, query_input] {
        if s != "stdin" && !is_pgi_input(s) {
            inputs.push(sibling_pgi_path(Path::new(s)).display().to_string());
        }
    }
    if let Some(s) = args.get_one::<String>("ref_seq") {
        inputs.push(s.clone());
    }
    if let Some(s) = args.get_one::<String>("query_seq") {
        inputs.push(s.clone());
    }
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, inputs.iter().map(|s| s.as_str()))?;
    let params = pgr::libs::pgi::align::AlignParams {
        freq: *args.get_one::<u32>("freq").unwrap(),
        min_shared: args.get_one::<usize>("min_shared").copied(),
    };
    let keep = args.get_flag("keep_index");
    let parallel = *args.get_one::<usize>("parallel").unwrap();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(parallel)
        .build()
        .context("building align thread pool")?;
    // Everything from the automatic index build through the merge/extension
    // runs inside the `--parallel` pool, so the option caps rayon usage for
    // the whole command (the sibling-index build used to run on the global
    // pool regardless of --parallel).
    let psls = pool.install(|| -> anyhow::Result<Vec<pgr::libs::fmt::psl::Psl>> {
        let mut tmp: Option<tempfile::TempDir> = None;
        let SideInput {
            index: ref_index,
            seqs: ref_side_seqs,
        } = resolve_side(args, ref_input, "reference", &mut tmp, keep)?;
        let SideInput {
            index: query_index,
            seqs: query_side_seqs,
        } = if self_mode {
            // Self-alignment resolves the same input once and reuses it on
            // both sides (the sequence copy is bounded by the input size).
            SideInput {
                index: ref_index.clone(),
                seqs: ref_side_seqs.clone(),
            }
        } else {
            resolve_side(args, query_input, "query", &mut tmp, keep)?
        };
        let _tmp_guard = tmp;

        // The reference index is consumed as a stream by the merge; the query
        // index is memory-mapped (FastGA's GIX model) and decoded on demand,
        // so neither index is materialized in full.
        let mut r1 = pgr::reader(&ref_index)?;
        let mut a = pgr::libs::pgi::PgiStream::open(&mut r1)?;
        let b = pgr::libs::pgi::PgiMmap::open(std::path::Path::new(&query_index))?;

        // Extension sequences come from genome inputs directly (validated
        // against the index, which matters when a sibling index was reused)
        // or from --ref-seq/--query-seq for .pgi inputs (validated the same
        // way).
        let mut ref_seqs = resolve_seqs(
            args,
            ref_side_seqs,
            a.header().contigs.as_slice(),
            "reference",
            "ref_seq",
        )?;
        let mut query_seqs =
            resolve_seqs(args, query_side_seqs, b.contigs(), "query", "query_seq")?;
        // Self mode aligns one input to itself, so a single --ref-seq (or
        // --query-seq) on a .pgi input supplies the extension sequences for
        // both sides (previously it errored unless both flags were given).
        if self_mode {
            if query_seqs.is_empty() && !ref_seqs.is_empty() {
                query_seqs = ref_seqs.clone();
            } else if ref_seqs.is_empty() && !query_seqs.is_empty() {
                ref_seqs = query_seqs.clone();
            }
        };
        if ref_seqs.is_empty() != query_seqs.is_empty() {
            anyhow::bail!(
                "extension sequences are needed for both sides (genome inputs, or \
                 --ref-seq/--query-seq for .pgi inputs)"
            );
        }
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
        "wrote {} PSL blocks (freq={}) to {}",
        psls.len(),
        params.freq,
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
    if cached.exists() && !sibling_index_stale(input, &cached) {
        let (ck, cs, cw) = read_index_params(&cached)?;
        let k = *args.get_one::<usize>("kmer").unwrap();
        let smer = *args.get_one::<usize>("smer").unwrap();
        let window = *args.get_one::<usize>("window").unwrap();
        // The current parameters (explicit or the defaults) must match the
        // cached index; a default `-k 40` run silently reusing a `k=20`
        // sibling index would report k=40 semantics with k=20 seeds.
        if k != ck {
            anyhow::bail!(
                "--kmer {k} conflicts with the cached index {} (k={ck})",
                cached.display()
            );
        }
        if smer != cs {
            anyhow::bail!(
                "--smer {smer} conflicts with the cached index {} (smer={cs})",
                cached.display()
            );
        }
        if window != cw {
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
    let idx = pgr::libs::pgi::build::build_from_seqs(seqs.clone(), k, smer, window, false, true)?;
    let out = if keep {
        cached.clone()
    } else {
        if tmp.is_none() {
            *tmp = Some(tempfile::TempDir::new().context("creating temporary index directory")?);
        }
        tmp.as_ref()
            .expect("temporary index directory initialized above")
            .path()
            .join(format!("{label}.pgi"))
    };
    // `idx.write` issues one write per occurrence record (millions); the
    // buffered writer turns those into large chunks instead of one syscall
    // per record (the standalone `pgr pgi build` path uses `pgr::writer`).
    let mut w = std::io::BufWriter::new(std::fs::File::create(&out)?);
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

/// Whether the genome input was modified after the sibling index was built.
///
/// The index stores only k-mers, so a FASTA edited in place (same contig
/// names/lengths, different sequence) would otherwise silently reuse the
/// stale index; rebuild instead (same mtime convention as the e-kmer cache).
fn sibling_index_stale(input: &str, index: &Path) -> bool {
    let (Ok(input_m), Ok(index_m)) = (
        std::fs::metadata(input).and_then(|m| m.modified()),
        std::fs::metadata(index).and_then(|m| m.modified()),
    ) else {
        return false;
    };
    input_m > index_m
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

/// Sibling index path for a genome input: ref.fa and ref.2bit map to ref.pgi;
/// a gzipped input (ref.fa.gz) gets its own ref.fa.pgi so it cannot silently
/// reuse the index of a plain ref.fa when the two files hold different
/// sequences (same contig names/lengths, different content).
fn sibling_pgi_path(path: &Path) -> PathBuf {
    let mut p = path.to_path_buf();
    if p.extension().and_then(|e| e.to_str()) == Some("gz") {
        // `ref.fa.gz` -> `ref.fa.pgi`: drop the `.gz` and keep the `.fa` so
        // the gzipped input has its own index, distinct from a plain
        // `ref.fa` (which maps to `ref.pgi`). A shared sibling would reuse
        // the wrong index when both files exist with the same contig
        // names/lengths but different sequences.
        p.set_extension("");
        return PathBuf::from(format!("{}.pgi", p.display()));
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
        // Read soft-mask blocks lowercase so the automatic-index masking
        // (build_from_seqs `mask=true`, FastGA `-M`) skips them exactly like
        // FASTA lowercase regions. Reading unmasked would index and align
        // masked repeats that a FASTA sibling silently skips.
        pgr::libs::pgi::build::read_2bit(path, true)
    } else {
        pgr::libs::pgi::build::read_fasta(path)
    }
}
