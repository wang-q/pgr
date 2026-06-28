use clap::*;
use indexmap::IndexMap;
use noodles_core::Position;
use noodles_fasta as fasta;
use pgr::libs::loc;
use pgr::libs::paf::cigar::CigarOp;
use pgr::libs::paf::index::{PafIndex, QueryResult};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::BufRead;
use std::num::NonZeroUsize;

use super::query;

// Load TSV mapping genome_name -> bgzf_fasta_path.
// Lines starting with '#' are comments; blank lines are skipped.
pub(crate) fn load_fasta_tsv(path: &str) -> anyhow::Result<IndexMap<String, String>> {
    let f = fs::File::open(path)?;
    let mut map = IndexMap::new();
    for line in std::io::BufReader::new(f).lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        anyhow::ensure!(
            fields.len() >= 2,
            "invalid TSV line '{line}': expected 2 tab-separated columns"
        );
        let name = fields[0].to_string();
        let fasta_path = fields[1].to_string();
        if map.insert(name.clone(), fasta_path).is_some() {
            anyhow::bail!("duplicate genome name in TSV: {name}");
        }
    }
    Ok(map)
}

// One opened BGZF FASTA file with its .loc index and a per-name record cache.
pub(crate) struct FastaEntry {
    reader: loc::Input,
    loc_of: IndexMap<String, (u64, usize)>,
    cache: lru::LruCache<String, fasta::Record>,
}

// Manages multiple BGZF FASTA files keyed by file path, with a name -> file
// mapping so multiple genome names can share one file (multi-chrom).
pub(crate) struct FastaStore {
    files: HashMap<String, FastaEntry>,
    name_to_file: HashMap<String, String>,
}

impl FastaStore {
    pub(crate) fn new(seq_to_file: &IndexMap<String, String>) -> anyhow::Result<Self> {
        let mut files = HashMap::new();
        let name_to_file: HashMap<String, String> = seq_to_file
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();

        // Unique file paths
        let unique_paths: HashSet<&String> = seq_to_file.values().collect();
        for path in unique_paths {
            let loc_file = format!("{path}.loc");
            if !std::path::Path::new(&loc_file).is_file() {
                loc::create_loc(path, &loc_file, true)?;
            }
            let loc_of = loc::load_loc(&loc_file)?;
            let reader = loc::Input::Bgzf(
                noodles_bgzf::io::indexed_reader::Builder::default().build_from_path(path)?,
            );
            let cache = lru::LruCache::new(NonZeroUsize::new(8).unwrap());
            files.insert(
                path.clone(),
                FastaEntry {
                    reader,
                    loc_of,
                    cache,
                },
            );
        }

        Ok(Self {
            files,
            name_to_file,
        })
    }

    // Fetch sequence [start, end) (0-based, half-open) and the total sequence
    // length. Caches the underlying FASTA record keyed by `name`.
    pub(crate) fn fetch_range(
        &mut self,
        name: &str,
        start: i32,
        end: i32,
    ) -> anyhow::Result<(Vec<u8>, usize)> {
        let path = self
            .name_to_file
            .get(name)
            .ok_or_else(|| anyhow::anyhow!("sequence '{name}' not in FASTA store"))?;
        let entry = self
            .files
            .get_mut(path)
            .ok_or_else(|| anyhow::anyhow!("file '{path}' not opened"))?;

        if !entry.cache.contains(name) {
            let record = loc::fetch_record(&mut entry.reader, &entry.loc_of, name)?;
            entry.cache.put(name.to_string(), record);
        }
        let record = entry.cache.get(name).unwrap();
        let total_len = record.sequence().len();

        // noodles Position is 1-based inclusive; our coords are 0-based half-open.
        let start_pos = Position::new(start as usize + 1).unwrap();
        let end_pos = Position::new(end as usize).unwrap();
        let slice = record
            .sequence()
            .slice(start_pos..=end_pos)
            .ok_or_else(|| anyhow::anyhow!("slice [{start},{end}) out of range for '{name}'"))?;

        Ok((slice.as_ref().to_vec(), total_len))
    }
}

// Reverse-complement a DNA byte slice (ACGTN-aware, case-preserving).
// Non-ACGTN bytes are passed through unchanged.
pub(crate) fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    fn comp(b: u8) -> u8 {
        match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            other => other,
        }
    }
    seq.iter().rev().map(|&b| comp(b)).collect()
}

// Build aligned strings (query, target) by walking CIGAR over [ts, te).
// `q_seq` covers query[qs..qe), `t_seq` covers target[ts..te).
// CIGAR origin is (rec_ts, rec_qs). Ops before `ts` are skipped (with partial
// skip for =/X/M/D); ops at/after `te` are stopped.
#[allow(clippy::too_many_arguments)]
fn build_maf_block(
    cigar: &[CigarOp],
    rec_ts: i32,
    rec_qs: i32,
    ts: i32,
    te: i32,
    qs: i32,
    q_seq: &[u8],
    t_seq: &[u8],
) -> (String, String) {
    let mut ct = rec_ts;
    let mut cq = rec_qs;
    let mut q_aln = String::new();
    let mut t_aln = String::new();

    for op in cigar {
        if ct >= te {
            break;
        }
        let td = op.target_delta() as i32;
        let qd = op.query_delta() as i32;
        let len = op.len() as i32;
        let next_ct = ct + td;

        match op.op() {
            '=' | 'X' | 'M' => {
                // Consume both query and target.
                let os = ct.max(ts);
                let oe = next_ct.min(te);
                if os < oe {
                    let skip_t = os - ct;
                    let take = oe - os;
                    let q_idx = (cq + skip_t - qs) as usize;
                    let t_idx = (os - ts) as usize;
                    for j in 0..take {
                        q_aln.push(q_seq[q_idx + j as usize] as char);
                        t_aln.push(t_seq[t_idx + j as usize] as char);
                    }
                }
                ct = next_ct;
                cq += len;
            }
            'I' => {
                // Consume query only (td == 0). Include if ct is within [ts, te).
                if ct >= ts && ct < te {
                    let q_idx = (cq - qs) as usize;
                    for j in 0..len {
                        q_aln.push(q_seq[q_idx + j as usize] as char);
                        t_aln.push('-');
                    }
                }
                cq += len;
            }
            'D' => {
                // Consume target only (qd == 0).
                let os = ct.max(ts);
                let oe = next_ct.min(te);
                if os < oe {
                    let t_idx = (os - ts) as usize;
                    let take = oe - os;
                    for j in 0..take {
                        q_aln.push('-');
                        t_aln.push(t_seq[t_idx + j as usize] as char);
                    }
                }
                ct = next_ct;
            }
            _ => {}
        }
        // suppress unused-assignment warning when qd == 0 (D op)
        let _ = qd;
    }

    (q_aln, t_aln)
}

// Output pairwise MAF blocks. Each QueryResult becomes one `a` block with two
// `s` lines (target first, query second). Sequences are fetched on demand via
// FastaStore; CIGAR is walked directly (no POA refinement).
#[allow(clippy::type_complexity)]
fn output_maf(
    idx: &PafIndex,
    all_results: &[((String, i32, i32), Vec<QueryResult>)],
    fasta_store: &mut FastaStore,
) -> anyhow::Result<()> {
    println!("##maf version=1");
    for (_, results) in all_results {
        for (query_id, q_iv, t_iv, cigar, rec_ts, rec_qs, strand) in results {
            let qname = idx.id_to_name(*query_id).unwrap_or("?");
            let tname = idx.id_to_name(t_iv.metadata).unwrap_or("?");

            let (qs, qe) = if q_iv.first <= q_iv.last {
                (q_iv.first, q_iv.last)
            } else {
                (q_iv.last, q_iv.first)
            };
            let (ts, te) = if t_iv.first <= t_iv.last {
                (t_iv.first, t_iv.last)
            } else {
                (t_iv.last, t_iv.first)
            };

            let (q_seq_fwd, q_src_size) = fasta_store.fetch_range(qname, qs, qe)?;
            let (t_seq, t_src_size) = fasta_store.fetch_range(tname, ts, te)?;

            // For '-' strand records: PAF query coords are on the forward
            // strand, but CIGAR describes alignment columns against the
            // reverse-complemented query. RC the fetched forward sequence
            // and walk CIGAR from offset 0 so column order matches.
            //
            // `q_seq_for_aln` is RC(forward[qs..qe)) and covers RC offset
            // [rec_qe - qe, rec_qe - qs) where rec_qe = rec_qs + aligned_q_len.
            // `build_maf_block` indexes q_seq via `(cq + skip_t) - qs_eff`, so
            // qs_eff must be the RC offset of the sub-interval start
            // (rec_qe - qe), not 0 — otherwise sub-interval queries index
            // past the end of q_seq. Full-overlap queries have qs_eff = 0.
            //
            // MAF `start` for '-' strand = srcSize - qe (position on forward
            // strand of the first displayed base, per MAF spec).
            let (q_seq_for_aln, rec_qs_eff, qs_eff, q_strand, q_start_maf) = if *strand == '-' {
                let rc = reverse_complement(&q_seq_fwd);
                let aligned_q_len: i32 = cigar.iter().map(|op| op.query_delta() as i32).sum();
                let rec_qe = *rec_qs + aligned_q_len;
                let rc_sub_start = rec_qe - qe;
                (rc, 0, rc_sub_start, '-', q_src_size as i32 - qe)
            } else {
                (q_seq_fwd, *rec_qs, qs, '+', qs)
            };

            let (q_aln, t_aln) = build_maf_block(
                cigar,
                *rec_ts,
                rec_qs_eff,
                ts,
                te,
                qs_eff,
                &q_seq_for_aln,
                &t_seq,
            );

            // size = number of non-gap bases
            let q_size = q_aln.chars().filter(|c| *c != '-').count();
            let t_size = t_aln.chars().filter(|c| *c != '-').count();

            println!("a");
            println!("s\t{tname}\t{ts}\t{t_size}\t+\t{t_src_size}\t{t_aln}");
            println!("s\t{qname}\t{q_start_maf}\t{q_size}\t{q_strand}\t{q_src_size}\t{q_aln}");
            println!();
        }
    }
    Ok(())
}

// One entry to feed into POA: aligned sequence plus metadata for the MAF `s` line.
pub(crate) struct MsaEntry {
    pub name: String,
    pub start: i32,      // MAF start (forward-strand coordinate)
    pub strand: char,    // '+' or '-'
    pub src_size: usize, // total sequence length
    pub seq: Vec<u8>,    // sequence in alignment orientation (already RC if '-')
}

// Collect target + query sequences for one region into MsaEntry list.
// Target is taken from the first result's t_iv; queries are RC'd if '-' strand.
// Skips a query that duplicates the target (BFS self-loop via mirror index).
pub(crate) fn build_msa_entries(
    idx: &PafIndex,
    tname_region: &str,
    results: &[QueryResult],
    fasta_store: &mut FastaStore,
) -> anyhow::Result<Vec<MsaEntry>> {
    let mut entries: Vec<MsaEntry> = Vec::with_capacity(results.len() + 1);

    // Target entry from the first result.
    let (_, _, t_iv_first, _, _, _, _) = &results[0];
    let tname = idx.id_to_name(t_iv_first.metadata).unwrap_or(tname_region);
    let (ts, te) = if t_iv_first.first <= t_iv_first.last {
        (t_iv_first.first, t_iv_first.last)
    } else {
        (t_iv_first.last, t_iv_first.first)
    };
    let (t_seq, t_src_size) = fasta_store.fetch_range(tname, ts, te)?;
    entries.push(MsaEntry {
        name: tname.to_string(),
        start: ts,
        strand: '+',
        src_size: t_src_size,
        seq: t_seq,
    });

    // Query entries. Skip a query that duplicates the target entry.
    let t_key = (tname.to_string(), ts, '+', t_src_size);
    for (query_id, q_iv, _t_iv, _cigar, _rec_ts, _rec_qs, strand) in results {
        let qname = idx.id_to_name(*query_id).unwrap_or("?");
        let (qs, qe) = if q_iv.first <= q_iv.last {
            (q_iv.first, q_iv.last)
        } else {
            (q_iv.last, q_iv.first)
        };
        let (q_seq_fwd, q_src_size) = fasta_store.fetch_range(qname, qs, qe)?;
        let (seq, start, strand_char) = if *strand == '-' {
            (reverse_complement(&q_seq_fwd), q_src_size as i32 - qe, '-')
        } else {
            (q_seq_fwd, qs, '+')
        };
        let q_key = (qname.to_string(), start, strand_char, q_src_size);
        if q_key == t_key {
            continue;
        }
        entries.push(MsaEntry {
            name: qname.to_string(),
            start,
            strand: strand_char,
            src_size: q_src_size,
            seq,
        });
    }
    Ok(entries)
}

// Output multi-way MAF blocks via POA. For each region, collect target +
// all query sequences (queries RC'd if '-' strand), feed them into the POA
// engine, and emit one `a` block with N `s` lines. CIGAR is ignored.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn output_maf_msa(
    idx: &PafIndex,
    all_results: &[((String, i32, i32), Vec<QueryResult>)],
    fasta_store: &mut FastaStore,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> anyhow::Result<()> {
    let params = pgr::libs::poa::AlignmentParams {
        match_score,
        mismatch_score,
        gap_open,
        gap_extend,
    };

    println!("##maf version=1");
    for ((tname_region, _, _), results) in all_results {
        if results.is_empty() {
            continue;
        }

        let entries = build_msa_entries(idx, tname_region, results, fasta_store)?;

        // Run POA MSA.
        let mut poa =
            pgr::libs::poa::Poa::new(params.clone(), pgr::libs::poa::AlignmentType::Global);
        for e in &entries {
            poa.add_sequence(&e.seq);
        }
        let msa = poa.msa();

        // Emit the MAF block.
        println!("a");
        for (e, aln) in entries.iter().zip(msa.iter()) {
            let size = aln.chars().filter(|c| *c != '-').count();
            println!(
                "s\t{}\t{}\t{}\t{}\t{}\t{}",
                e.name, e.start, size, e.strand, e.src_size, aln
            );
        }
        println!();
    }
    Ok(())
}

pub fn make_subcommand() -> Command {
    query::add_query_args(
        Command::new("to-maf")
            .arg(
                Arg::new("fasta_tsv")
                    .long("fasta-tsv")
                    .short('f')
                    .required(true)
                    .num_args(1)
                    .help("TSV file: genome_name <tab> bgzf_fasta_path"),
            )
            .arg(
                Arg::new("msa")
                    .long("msa")
                    .num_args(0)
                    .help("Merge results per region into a multi-way MAF block via POA"),
            )
            .arg(
                Arg::new("match_score")
                    .long("match")
                    .num_args(1)
                    .default_value("5")
                    .value_parser(clap::value_parser!(i32))
                    .allow_negative_numbers(true)
                    .help("POA match score (default: 5)"),
            )
            .arg(
                Arg::new("mismatch_score")
                    .long("mismatch")
                    .num_args(1)
                    .default_value("-4")
                    .value_parser(clap::value_parser!(i32))
                    .allow_negative_numbers(true)
                    .help("POA mismatch score (default: -4)"),
            )
            .arg(
                Arg::new("gap_open")
                    .long("gap-open")
                    .num_args(1)
                    .default_value("-8")
                    .value_parser(clap::value_parser!(i32))
                    .allow_negative_numbers(true)
                    .help("POA gap open penalty (default: -8)"),
            )
            .arg(
                Arg::new("gap_extend")
                    .long("gap-extend")
                    .num_args(1)
                    .default_value("-6")
                    .value_parser(clap::value_parser!(i32))
                    .allow_negative_numbers(true)
                    .help("POA gap extend penalty (default: -6)"),
            ),
    )
    .about("Query PAF index and output pairwise or multi-way MAF")
    .after_help(
        r###"
Queries a PAF file or saved index (same logic as `pgr paf query`) and
outputs MAF blocks.

Default mode (pairwise): each query result becomes one 2-sequence MAF
block restored directly from CIGAR. Alignments are assumed to be
already refined by chain/net — no POA refinement is performed.

--msa mode (multi-way): merge all query results of each region into a
single multi-sequence MAF block via POA. Sequences (target first, then
each query, '-' strand reverse-complemented) are fed into the POA
engine; CIGAR is ignored. Best used with --transitive to gather all
homologous fragments of a region.

-f/--fasta-tsv (required):
  TSV with two columns: genome_name <tab> bgzf_fasta_path
  Each genome_name must match a query/target name in the PAF index.
  A FASTA file may be referenced by multiple genome_names (multi-chrom).
  All genome names in the PAF index must be present in the TSV (strict
  validation — missing entries cause a hard error).

Notes:
* Input PAF files should contain cg:Z: tags for accurate projection
* Supports both plain text and gzipped (.gz) files (including BGZF)
* Reads from stdin if input file is 'stdin'

Examples:
1. Single region to pairwise MAF:
   pgr paf to-maf alignments.paf chr1:1000-5000 -f genomes.tsv

2. Multi-way MSA with transitive BFS:
   pgr paf to-maf alignments.paf chr1:1000-5000 -t --msa -f genomes.tsv

3. Batch query from BED regions:
   pgr paf to-maf alignments.paf.idx -b regions.bed -f genomes.tsv

"###,
    )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let tsv_path = args.get_one::<String>("fasta_tsv").unwrap();
    let seq_to_file = load_fasta_tsv(tsv_path)?;

    let (idx, all_results) = query::run_query(args)?;

    // Validate: every name in the PAF index must be present in the TSV.
    let mut missing: Vec<&str> = idx
        .names
        .keys()
        .filter(|n| !seq_to_file.contains_key(*n))
        .map(|n| n.as_str())
        .collect();
    missing.sort_unstable();
    if !missing.is_empty() {
        anyhow::bail!(
            "FASTA TSV is missing {} genome(s) present in PAF index: {}",
            missing.len(),
            missing.join(", ")
        );
    }

    let mut fasta_store = FastaStore::new(&seq_to_file)?;
    if args.get_flag("msa") {
        let match_score = *args.get_one::<i32>("match_score").unwrap();
        let mismatch_score = *args.get_one::<i32>("mismatch_score").unwrap();
        let gap_open = *args.get_one::<i32>("gap_open").unwrap();
        let gap_extend = *args.get_one::<i32>("gap_extend").unwrap();
        output_maf_msa(
            &idx,
            &all_results,
            &mut fasta_store,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
        )?;
    } else {
        output_maf(&idx, &all_results, &mut fasta_store)?;
    }
    Ok(())
}
