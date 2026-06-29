use crate::libs::paf::cigar::CigarOp;
use crate::libs::paf::fasta::FastaStore;
use crate::libs::paf::index::{PafIndex, QueryResult};

/// Reverse-complement a DNA byte slice (ACGTN-aware, case-preserving).
/// Non-ACGTN bytes are passed through unchanged.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
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
