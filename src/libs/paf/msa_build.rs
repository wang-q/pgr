use crate::libs::alignment::coords;
use crate::libs::nt;
use crate::libs::paf::cigar::CigarOp;
use crate::libs::paf::fasta::FastaStore;
use crate::libs::paf::index::{PafIndex, QueryResult};
use crate::libs::poa::{self, AlignmentParams};
use std::collections::HashSet;

/// Return `(start, end)` with `start <= end` from an oriented interval.
/// PAF intervals may be stored as `(first, last)` in either order; this
/// normalizes them to ascending half-open ranges used by all emitters.
pub fn orient_interval(first: i32, last: i32) -> (i32, i32) {
    if first <= last {
        (first, last)
    } else {
        (last, first)
    }
}

/// Build aligned strings (query, target) by walking CIGAR over [ts, te).
/// `q_seq` covers query[qs..qe), `t_seq` covers target[ts..te).
/// CIGAR origin is (rec_ts, rec_qs). Ops before `ts` are skipped (with partial
/// skip for =/X/M/D); ops at/after `te` are stopped.
///
/// Returns `Err` if the CIGAR is inconsistent with the provided sequence
/// lengths (defensive guard against malformed index data).
#[allow(clippy::too_many_arguments)]
pub fn build_maf_block(
    cigar: &[CigarOp],
    rec_ts: i32,
    rec_qs: i32,
    ts: i32,
    te: i32,
    qs: i32,
    q_seq: &[u8],
    t_seq: &[u8],
) -> anyhow::Result<(String, String)> {
    let mut ct = rec_ts;
    let mut cq = rec_qs;
    let mut q_aln = String::new();
    let mut t_aln = String::new();
    let q_len = q_seq.len() as i32;
    let t_len = t_seq.len() as i32;

    for op in cigar {
        if ct >= te {
            break;
        }
        let td = op.target_delta() as i32;
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
                    let q_idx = cq + skip_t - qs;
                    let t_idx = os - ts;
                    if q_idx < 0 || t_idx < 0 || q_idx + take > q_len || t_idx + take > t_len {
                        anyhow::bail!(
                            "build_maf_block: CIGAR/sequence mismatch at M op (cq={}, skip_t={}, qs={}, q_idx={}, take={}, q_len={}, t_idx={}, t_len={})",
                            cq, skip_t, qs, q_idx, take, q_len, t_idx, t_len
                        );
                    }
                    let q_idx = q_idx as usize;
                    let t_idx = t_idx as usize;
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
                    let q_idx = cq - qs;
                    if q_idx < 0 || q_idx + len > q_len {
                        anyhow::bail!(
                            "build_maf_block: CIGAR/sequence mismatch at I op (cq={}, qs={}, q_idx={}, len={}, q_len={})",
                            cq, qs, q_idx, len, q_len
                        );
                    }
                    let q_idx = q_idx as usize;
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
                    let t_idx = os - ts;
                    let take = oe - os;
                    if t_idx < 0 || t_idx + take > t_len {
                        anyhow::bail!(
                            "build_maf_block: CIGAR/sequence mismatch at D op (os={}, ts={}, t_idx={}, take={}, t_len={})",
                            os, ts, t_idx, take, t_len
                        );
                    }
                    let t_idx = t_idx as usize;
                    for j in 0..take {
                        q_aln.push('-');
                        t_aln.push(t_seq[t_idx + j as usize] as char);
                    }
                }
                ct = next_ct;
            }
            _ => {}
        }
    }

    Ok((q_aln, t_aln))
}

/// One entry to feed into POA: aligned sequence plus metadata for the MAF `s` line.
pub struct MsaEntry {
    pub name: String,
    pub start: i32,      // MAF start (forward-strand coordinate)
    pub strand: char,    // '+' or '-'
    pub src_size: usize, // total sequence length
    pub seq: Vec<u8>,    // sequence in alignment orientation (already RC if '-')
}

/// Collect target + query sequences for one region into MsaEntry list.
/// Target is taken from the first result's t_iv; queries are RC'd if '-' strand.
/// Skips a query that duplicates the target (BFS self-loop via mirror index).
pub fn build_msa_entries(
    idx: &PafIndex,
    tname_region: &str,
    results: &[QueryResult],
    fasta_store: &mut FastaStore,
) -> anyhow::Result<Vec<MsaEntry>> {
    let mut entries: Vec<MsaEntry> = Vec::with_capacity(results.len() + 1);

    // Target entry from the first result.
    let (_, _, t_iv_first, _, _, _, _) = &results[0];
    let tname = idx.id_to_name(t_iv_first.metadata).unwrap_or(tname_region);
    let (ts, te) = orient_interval(t_iv_first.first, t_iv_first.last);
    let (t_seq, t_src_size) = fasta_store.fetch_range(tname, ts, te)?;
    entries.push(MsaEntry {
        name: tname.to_string(),
        start: ts,
        strand: '+',
        src_size: t_src_size,
        seq: t_seq,
    });

    // Transitive BFS can report the same query fragment once per alignment
    // path (10 genomes -> dozens of near-identical copies). Deduplicate exact
    // duplicates, then merge overlapping same-strand intervals of the same
    // query, so POA receives one entry per fragment instead of N copies.
    let t_key = (tname.to_string(), ts, '+', t_src_size);
    let mut seen: HashSet<(String, i32, i32, char)> = HashSet::new();
    let mut intervals: Vec<(String, i32, i32, char)> = Vec::new();
    for (query_id, q_iv, _t_iv, _cigar, _rec_ts, _rec_qs, strand) in results {
        let qname = idx.id_to_name(*query_id).unwrap_or("?");
        let (qs, qe) = orient_interval(q_iv.first, q_iv.last);
        if seen.insert((qname.to_string(), qs, qe, *strand)) {
            intervals.push((qname.to_string(), qs, qe, *strand));
        }
    }
    // Group by (name, strand) so same-strand intervals are contiguous for the
    // overlap merge below; a '-' strand copy of the same region would
    // otherwise sit between two '+' fragments and block their merge.
    intervals.sort_by(|a, b| a.0.cmp(&b.0).then(a.3.cmp(&b.3)).then(a.1.cmp(&b.1)));
    let mut merged: Vec<(String, i32, i32, char)> = Vec::new();
    for (qname, qs, qe, strand) in intervals {
        if let Some(last) = merged.last_mut() {
            // Merge overlapping OR exactly-touching same-strand intervals.
            // Exact-touch (qs == last.2) happens when a target deletion splits
            // one continuous query region into two abutting fragments; keeping
            // them separate would let the per-name dedup below silently drop
            // the second fragment's bases.
            if last.0 == qname && last.3 == strand && qs <= last.2 {
                last.2 = last.2.max(qe);
                continue;
            }
        }
        merged.push((qname, qs, qe, strand));
    }

    // Query entries. One block per region -> each genome appears at most once:
    // MAF blocks require unique sequence names and VCF sample columns must not
    // repeat. Genome-internal duplicate loci (e.g. rrn operons, origin-crossing
    // wraparound) surface as extra entries; keep the first per name (sorted
    // '+' before '-' for the same name, so the near-reference-strand copy wins).
    // Skip a query that duplicates the target entry.
    let mut seen_names: HashSet<String> = HashSet::new();
    for (qname, qs, qe, strand) in merged {
        if !seen_names.insert(qname.clone()) {
            continue;
        }
        let (q_seq_fwd, q_src_size) = fasta_store.fetch_range(&qname, qs, qe)?;
        let (seq, start, strand_char) = if strand == '-' {
            (
                nt::rev_comp(&q_seq_fwd).collect::<Vec<u8>>(),
                coords::reverse_range_pair(qs, qe, q_src_size as i32).0,
                '-',
            )
        } else {
            (q_seq_fwd, qs, '+')
        };
        let q_key = (qname.clone(), start, strand_char, q_src_size);
        if q_key == t_key {
            continue;
        }
        entries.push(MsaEntry {
            name: qname,
            start,
            strand: strand_char,
            src_size: q_src_size,
            seq,
        });
    }
    Ok(entries)
}

/// One pairwise alignment record restored from a CIGAR.
/// Carries aligned strings plus the metadata needed by MAF / FAS emitters.
pub struct PairwiseBlock {
    pub qname: String,
    pub tname: String,
    pub q_aln: String,
    pub t_aln: String,
    /// MAF `start` field: forward-strand coord of first displayed base.
    /// '+' strand: qs. '-' strand: src_size - q_end_fwd.
    pub q_start_maf: i32,
    pub q_strand: char,
    pub q_src_size: usize,
    pub t_start: i32,
    pub t_end: i32,
    pub t_src_size: usize,
}

/// Project one `QueryResult` through its CIGAR and fetch sequences,
/// returning alignment strings plus metadata for both MAF and FAS emitters.
///
/// For `-` strand records: PAF query coords are on the forward strand, but
/// CIGAR describes alignment columns against the reverse-complemented query.
/// We RC the fetched forward sequence and walk CIGAR from offset 0 so column
/// order matches. `q_start_maf` is set to `src_size - qe` per MAF spec.
pub fn build_pairwise_block(
    idx: &PafIndex,
    result: &QueryResult,
    fasta_store: &mut FastaStore,
) -> anyhow::Result<PairwiseBlock> {
    let (query_id, q_iv, t_iv, cigar, rec_ts, rec_qs, strand) = result;
    let qname = idx.id_to_name(*query_id).unwrap_or("?").to_string();
    let tname = idx.id_to_name(t_iv.metadata).unwrap_or("?").to_string();

    let (qs, qe) = orient_interval(q_iv.first, q_iv.last);
    let (ts, te) = orient_interval(t_iv.first, t_iv.last);

    let (q_seq_fwd, q_src_size) = fasta_store.fetch_range(&qname, qs, qe)?;
    let (t_seq, t_src_size) = fasta_store.fetch_range(&tname, ts, te)?;

    let (q_seq_for_aln, rec_qs_eff, qs_eff, q_strand, q_start_maf) = if *strand == '-' {
        let rc = nt::rev_comp(&q_seq_fwd).collect::<Vec<u8>>();
        let aligned_q_len: i32 = cigar.iter().map(|op| op.query_delta() as i32).sum();
        let rec_qe = *rec_qs + aligned_q_len;
        let rc_sub_start = rec_qe - qe;
        (
            rc,
            0,
            rc_sub_start,
            '-',
            coords::reverse_range_pair(qs, qe, q_src_size as i32).0,
        )
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
    )?;

    Ok(PairwiseBlock {
        qname,
        tname,
        q_aln,
        t_aln,
        q_start_maf,
        q_strand,
        q_src_size,
        t_start: ts,
        t_end: te,
        t_src_size,
    })
}

/// Run POA global MSA on a slice of `MsaEntry` and return one aligned string
/// per entry (parallel order). Thin wrapper around `Poa::new` + `add_sequence`
/// + `msa()` used by both `to-fas --msa` and `to-maf --msa`.
pub fn run_poa_msa(entries: &[MsaEntry], params: AlignmentParams) -> Vec<String> {
    let mut poa = poa::Poa::new(params, poa::AlignmentType::Global);
    for e in entries {
        poa.add_sequence(&e.seq);
    }
    poa.msa()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::paf::index::PafIndex;
    use coitrees::Interval;
    use indexmap::IndexMap;
    use std::collections::HashMap;
    use std::io::Write;

    /// Write a single-sequence BGZF FASTA and return its path (with .gzi built).
    fn write_bgzf_fasta(dir: &std::path::Path, name: &str, seq: &str) -> String {
        let path = dir.join(format!("{name}.fa.gz"));
        let file = std::fs::File::create(&path).unwrap();
        let mut writer = crate::libs::bgzf::BgzfWriter::new(file).unwrap();
        writeln!(writer, ">{name}").unwrap();
        writeln!(writer, "{seq}").unwrap();
        writer.finish().unwrap();
        crate::libs::fmt::fa::build_gzi_index(path.to_str().unwrap()).unwrap();
        path.to_string_lossy().to_string()
    }

    #[test]
    fn test_build_msa_entries_merges_touching_fragments() {
        let dir = tempfile::tempdir().unwrap();
        let q_path = write_bgzf_fasta(dir.path(), "Q", &"A".repeat(300));
        let t_path = write_bgzf_fasta(dir.path(), "T", &"C".repeat(300));
        let mut map = IndexMap::new();
        map.insert("Q".to_string(), q_path);
        map.insert("T".to_string(), t_path);
        let mut store = super::super::fasta::FastaStore::new(&map).unwrap();

        let mut names = IndexMap::new();
        names.insert("Q".to_string(), 0u32);
        names.insert("T".to_string(), 1u32);
        let idx = PafIndex {
            names,
            trees: HashMap::new(),
            reverse_trees: HashMap::new(),
            lazy_source: None,
            lazy_source_path: None,
        };

        // Two abutting same-strand fragments of Q (a target deletion splits one
        // continuous query region into [0,100) and [100,200)). They must merge
        // into a single [0,200) entry rather than being dropped by the
        // per-name dedup.
        let mk = |q_start: i32, q_end: i32| {
            (
                0u32, // query_id = Q
                Interval::new(q_start, q_end, 0u32),
                Interval::new(0, 100, 1u32), // target T [0,100)
                vec![CigarOp::new((q_end - q_start) as u32, '=')],
                0,
                0,
                '+',
            )
        };
        let results = vec![mk(0, 100), mk(100, 200)];

        let entries = build_msa_entries(&idx, "T", &results, &mut store).unwrap();
        // Target entry + ONE merged query entry covering [0,200).
        assert_eq!(
            entries.len(),
            2,
            "expected target + one merged query entry, got {}",
            entries.len()
        );
        let q_entry = entries.iter().find(|e| e.name == "Q").unwrap();
        assert_eq!(q_entry.start, 0);
        assert_eq!(
            q_entry.seq.len(),
            200,
            "merged query must span [0,200), got {} bases",
            q_entry.seq.len()
        );
    }
}
