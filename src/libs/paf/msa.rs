use crate::libs::nt;
use crate::libs::paf::cigar::CigarOp;
use crate::libs::paf::fasta::FastaStore;
use crate::libs::paf::index::{PafIndex, QueryResult};
use crate::libs::poa::{self, AlignmentParams};

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

    // Query entries. Skip a query that duplicates the target entry.
    let t_key = (tname.to_string(), ts, '+', t_src_size);
    for (query_id, q_iv, _t_iv, _cigar, _rec_ts, _rec_qs, strand) in results {
        let qname = idx.id_to_name(*query_id).unwrap_or("?");
        let (qs, qe) = orient_interval(q_iv.first, q_iv.last);
        let (q_seq_fwd, q_src_size) = fasta_store.fetch_range(qname, qs, qe)?;
        let (seq, start, strand_char) = if *strand == '-' {
            (
                nt::rev_comp(&q_seq_fwd).collect::<Vec<u8>>(),
                q_src_size as i32 - qe,
                '-',
            )
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

/// One pairwise alignment record restored from a CIGAR.
/// Carries aligned strings plus the metadata needed by MAF / FAS emitters.
pub struct PairwiseBlock {
    pub qname: String,
    pub tname: String,
    pub q_aln: String,
    pub t_aln: String,
    /// Forward-strand start (qs in PAF coords). Used by FAS emitter.
    pub q_start_fwd: i32,
    /// Forward-strand end (qe in PAF coords). Used by FAS emitter.
    pub q_end_fwd: i32,
    /// MAF `start` field: forward-strand coord of first displayed base.
    /// '+' strand: == q_start_fwd. '-' strand: src_size - q_end_fwd.
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

    Ok(PairwiseBlock {
        qname,
        tname,
        q_aln,
        t_aln,
        q_start_fwd: qs,
        q_end_fwd: qe,
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
