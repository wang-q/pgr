//! VCF output from POA multiple sequence alignment of PAF query results.

use crate::libs::paf::fasta::FastaStore;
use crate::libs::paf::index::{PafIndex, QueryResult};
use crate::libs::paf::msa::build_msa_entries;
use crate::libs::poa;
use std::io::Write;

// Emit one VCF row. ref_allele is REF; alt_alleles are distinct non-REF
// alleles (joined by ','); sample_alleles[i] is sample i's allele string
// (empty or non-ACGT -> GT='.'). GT: 0=REF, 1..=N=ALT index, '.'=other.
fn emit_vcf_row<W: Write>(
    writer: &mut W,
    chrom: &str,
    pos: i32,
    ref_allele: &str,
    alt_alleles: &[String],
    sample_alleles: &[String],
) -> anyhow::Result<()> {
    let alt = alt_alleles.join(",");
    let mut row = String::new();
    row.push_str(chrom);
    row.push('\t');
    row.push_str(&pos.to_string());
    row.push_str("\t.\t");
    row.push_str(ref_allele);
    row.push('\t');
    row.push_str(&alt);
    row.push_str("\t.\t.\t.\tGT");
    for allele in sample_alleles {
        row.push('\t');
        let gt = if allele.is_empty() || allele == "-" {
            ".".to_string()
        } else if allele == ref_allele {
            "0".to_string()
        } else {
            match alt_alleles.iter().position(|a| a == allele) {
                Some(i) => (i + 1).to_string(),
                None => ".".to_string(),
            }
        };
        row.push_str(&gt);
    }
    row.push('\n');
    writer.write_all(row.as_bytes())?;
    Ok(())
}

// Left-align a set of indels (INS or DEL) by shifting the anchor position
// leftward as far as the reference allows. All non-empty indel sequences
// shift together: a step is taken only when the reference base immediately
// before the anchor equals the last base of *every* non-empty indel seq.
//
// `anchor_offset` is the index of the anchor base within `target_ext`
// (the extended target = prefix + aligned region, all uppercase).
// `indel_seqs[i]` is sample i's inserted (INS) or deleted (DEL) sequence;
// empty entries are samples with no indel at this site (they stay empty).
//
// Returns (new_anchor_offset, new_anchor_base, new_indel_seqs) where each
// new_indel_seqs[i] is the left-aligned version of indel_seqs[i].
fn left_align_indels(
    anchor_offset: usize,
    anchor_base: u8,
    indel_seqs: &[Vec<u8>],
    target_ext: &[u8],
) -> (usize, u8, Vec<Vec<u8>>) {
    let mut offset = anchor_offset;
    let mut anchor = anchor_base;
    let mut seqs: Vec<Vec<u8>> = indel_seqs.to_vec();

    loop {
        if offset == 0 {
            break;
        }
        let prev_base = target_ext[offset - 1];
        // All non-empty seqs must have last base == prev_base to shift.
        let can_shift = seqs
            .iter()
            .filter(|s| !s.is_empty())
            .all(|s| s.last().map(|&b| b == prev_base).unwrap_or(false));
        if !can_shift {
            break;
        }
        // Shift: each non-empty seq pops its last base and pushes the old
        // anchor. The new anchor becomes prev_base.
        for s in seqs.iter_mut() {
            if !s.is_empty() {
                s.pop();
                s.push(anchor);
            }
        }
        anchor = prev_base;
        offset -= 1;
    }

    (offset, anchor, seqs)
}

// Output VCF records from POA MSA of each region. Three variant classes
// are emitted: substitutions (single target non-gap column with ≥1 differing
// query), INS (consecutive target gap columns; REF = 1bp anchor at the
// preceding non-gap column, ALT = anchor + inserted bases), and DEL
// (consecutive target non-gap columns where ≥1 query has gap; REF = anchor
// + target segment, ALT = anchor + per-query non-gap bases; a fully-deleted
// sample gets ALT = anchor). INS and DEL are left-aligned against the
// reference: the anchor position shifts leftward while the reference base
// before the anchor equals the last base of every non-empty indel seq.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn output_vcf<W: Write>(
    writer: &mut W,
    idx: &PafIndex,
    all_results: &[((String, i32, i32), Vec<QueryResult>)],
    fasta_store: &mut FastaStore,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> anyhow::Result<()> {
    let params = poa::AlignmentParams {
        match_score,
        mismatch_score,
        gap_open,
        gap_extend,
    };

    let mut header_written = false;

    for ((tname_region, _, _), results) in all_results {
        if results.is_empty() {
            continue;
        }

        let entries = build_msa_entries(idx, tname_region, results, fasta_store)?;

        // Run POA MSA.
        let mut poa = poa::Poa::new(params.clone(), poa::AlignmentType::Global);
        for e in &entries {
            poa.add_sequence(&e.seq);
        }
        let msa = poa.msa();

        if msa.is_empty() {
            continue;
        }

        // Write VCF header once, using the entry names as sample columns.
        if !header_written {
            writer.write_all(b"##fileformat=VCFv4.2\n")?;
            writer
                .write_all(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")?;
            let mut header = String::from("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
            for e in &entries {
                header.push('\t');
                header.push_str(&e.name);
            }
            header.push('\n');
            writer.write_all(header.as_bytes())?;
            header_written = true;
        }

        let n_seq = entries.len();
        let aln_len = msa[0].len();
        let target = &entries[0];

        // Build extended target = prefix + aligned region (all uppercase) for
        // indel left-alignment. prefix covers up to 1000bp before the region
        // so the anchor can be shifted left within the reference.
        const LEFT_ALIGN_BUFFER: i32 = 1000;
        let prefix_start = (target.start - LEFT_ALIGN_BUFFER).max(0);
        let (target_prefix, _) = if target.start > prefix_start {
            fasta_store.fetch_range(&target.name, prefix_start, target.start)?
        } else {
            (Vec::new(), 0)
        };
        let mut target_ext: Vec<u8> = target_prefix;
        target_ext.extend_from_slice(&target.seq);
        for b in target_ext.iter_mut() {
            *b = b.to_ascii_uppercase();
        }
        // t_start_offset: index of target.start within target_ext.
        let t_start_offset = (target.start - prefix_start) as usize;

        // Walk MSA columns with a while loop so we can advance past indel
        // regions. t_aln_pos counts target non-gap columns processed (used
        // to derive VCF POS from target.start). Three cases:
        //   1. INS: consecutive target gap columns (REF = 1bp anchor at the
        //      preceding non-gap column, ALT = anchor + inserted bases).
        //   2. DEL: consecutive target non-gap columns where ≥1 query has
        //      gap (REF = target segment, ALT = per-query non-gap concat).
        //   3. SNP: single target non-gap column with no gaps.
        let mut col: usize = 0;
        let mut t_aln_pos: i32 = 0;
        while col < aln_len {
            let t_base = msa[0].as_bytes()[col];
            if t_base == b'-' {
                // INS region: collect consecutive target gap columns.
                let col_start = col;
                while col < aln_len && msa[0].as_bytes()[col] == b'-' {
                    col += 1;
                }
                let col_end = col;
                // Anchor = previous non-gap target column. Skip if none.
                if col_start == 0 {
                    continue;
                }
                let anchor_byte = msa[0].as_bytes()[col_start - 1];
                if anchor_byte == b'-' {
                    continue;
                }
                // Per-sample inserted bases (gaps and non-ACGT dropped).
                let sample_inserted: Vec<Vec<u8>> = msa
                    .iter()
                    .take(n_seq)
                    .map(|seq| {
                        (col_start..col_end)
                            .filter_map(|c| {
                                let b = seq.as_bytes()[c].to_ascii_uppercase();
                                if matches!(b, b'A' | b'C' | b'G' | b'T') {
                                    Some(b)
                                } else {
                                    None
                                }
                            })
                            .collect()
                    })
                    .collect();
                if !sample_inserted.iter().any(|s| !s.is_empty()) {
                    continue;
                }
                // Left-align: anchor offset in target_ext = t_aln_pos - 1
                // + t_start_offset (anchor is the last target non-gap column
                // before the INS gap, i.e. the (t_aln_pos)-th, 0-based idx
                // t_aln_pos - 1 in the aligned region).
                let anchor_offset = (t_aln_pos - 1 + t_start_offset as i32) as usize;
                let anchor_base = target_ext[anchor_offset];
                let (new_offset, new_anchor, new_inserted) =
                    left_align_indels(anchor_offset, anchor_base, &sample_inserted, &target_ext);
                let ref_allele = String::from(new_anchor as char);
                let sample_alleles: Vec<String> = new_inserted
                    .iter()
                    .map(|ins| {
                        let mut s = String::from(new_anchor as char);
                        for &b in ins {
                            s.push(b as char);
                        }
                        s
                    })
                    .collect();
                let mut alt_alleles: Vec<String> = Vec::new();
                for a in &sample_alleles {
                    if a != &ref_allele && !alt_alleles.contains(a) {
                        alt_alleles.push(a.clone());
                    }
                }
                if alt_alleles.is_empty() {
                    continue;
                }
                // POS (1-based) = prefix_start + new_offset + 1.
                let pos = prefix_start + new_offset as i32 + 1;
                emit_vcf_row(
                    writer,
                    &target.name,
                    pos,
                    &ref_allele,
                    &alt_alleles,
                    &sample_alleles,
                )?;
            } else {
                // target non-gap: check if any query has a gap here.
                let col_has_gap = msa.iter().take(n_seq).any(|s| s.as_bytes()[col] == b'-');
                if col_has_gap {
                    // DEL region: collect consecutive target non-gap columns
                    // where ≥1 query has a gap.
                    let col_start = col;
                    while col < aln_len {
                        let tb = msa[0].as_bytes()[col];
                        if tb == b'-' {
                            break;
                        }
                        let cg = msa.iter().take(n_seq).any(|s| s.as_bytes()[col] == b'-');
                        if cg {
                            col += 1;
                        } else {
                            break;
                        }
                    }
                    let col_end = col;
                    // Anchor = previous non-gap target column. Skip if none
                    // (can't represent a deletion without a 1bp anchor in VCF).
                    if col_start == 0 {
                        t_aln_pos += (col_end - col_start) as i32;
                        continue;
                    }
                    let anchor_byte = msa[0].as_bytes()[col_start - 1];
                    if anchor_byte == b'-' {
                        t_aln_pos += (col_end - col_start) as i32;
                        continue;
                    }
                    // target_segment: target bases in [col_start, col_end).
                    let target_segment: Vec<u8> = (col_start..col_end)
                        .filter_map(|c| {
                            let b = msa[0].as_bytes()[c].to_ascii_uppercase();
                            if matches!(b, b'A' | b'C' | b'G' | b'T') {
                                Some(b)
                            } else {
                                None
                            }
                        })
                        .collect();
                    if target_segment.is_empty() {
                        t_aln_pos += (col_end - col_start) as i32;
                        continue;
                    }
                    // Left-align target_segment (sole non-empty indel seq).
                    let anchor_offset = (t_aln_pos - 1 + t_start_offset as i32) as usize;
                    let anchor_base = target_ext[anchor_offset];
                    let indel_seqs = vec![target_segment.clone()];
                    let (new_offset, new_anchor, new_segs) =
                        left_align_indels(anchor_offset, anchor_base, &indel_seqs, &target_ext);
                    let new_segment = &new_segs[0];
                    // REF = new_anchor + new_segment.
                    let mut ref_allele = String::from(new_anchor as char);
                    for &b in new_segment {
                        ref_allele.push(b as char);
                    }
                    // Per-sample allele. Fully-deleted -> anchor only;
                    // fully-present -> REF; partial deletion -> best-effort
                    // (left-aligned anchor + original non-gap bases).
                    let sample_alleles: Vec<String> = msa
                        .iter()
                        .take(n_seq)
                        .map(|seq| {
                            let all_gap = (col_start..col_end).all(|c| seq.as_bytes()[c] == b'-');
                            if all_gap {
                                String::from(new_anchor as char)
                            } else {
                                let mut s = String::from(new_anchor as char);
                                for c in col_start..col_end {
                                    let b = seq.as_bytes()[c].to_ascii_uppercase();
                                    if matches!(b, b'A' | b'C' | b'G' | b'T') {
                                        s.push(b as char);
                                    }
                                }
                                s
                            }
                        })
                        .collect();
                    // POS (1-based) = prefix_start + new_offset + 1.
                    let pos = prefix_start + new_offset as i32 + 1;
                    t_aln_pos += (col_end - col_start) as i32;
                    let mut alt_alleles: Vec<String> = Vec::new();
                    for a in &sample_alleles {
                        if a != &ref_allele && !alt_alleles.contains(a) {
                            alt_alleles.push(a.clone());
                        }
                    }
                    if alt_alleles.is_empty() {
                        continue;
                    }
                    emit_vcf_row(
                        writer,
                        &target.name,
                        pos,
                        &ref_allele,
                        &alt_alleles,
                        &sample_alleles,
                    )?;
                } else {
                    // SNP: single target non-gap column, no gaps.
                    let ref_base = t_base.to_ascii_uppercase();
                    let ref_allele = String::from(ref_base as char);
                    let sample_alleles: Vec<String> = msa
                        .iter()
                        .take(n_seq)
                        .map(|seq| {
                            let b = seq.as_bytes()[col].to_ascii_uppercase();
                            if matches!(b, b'A' | b'C' | b'G' | b'T') {
                                String::from(b as char)
                            } else {
                                String::new()
                            }
                        })
                        .collect();
                    let mut alt_alleles: Vec<String> = Vec::new();
                    for a in &sample_alleles {
                        if !a.is_empty() && a != &ref_allele && !alt_alleles.contains(a) {
                            alt_alleles.push(a.clone());
                        }
                    }
                    let pos = target.start + t_aln_pos + 1;
                    t_aln_pos += 1;
                    col += 1;
                    if alt_alleles.is_empty() {
                        continue;
                    }
                    emit_vcf_row(
                        writer,
                        &target.name,
                        pos,
                        &ref_allele,
                        &alt_alleles,
                        &sample_alleles,
                    )?;
                }
            }
        }
    }

    Ok(())
}
