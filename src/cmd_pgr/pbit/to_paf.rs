//! Export PAF alignments embedded in a pbit archive.

use anyhow::Context;
use clap::{ArgMatches, Command};
use pgr::libs::paf::cigar::{block_identity, gap_compressed_identity, CigarOp};
use pgr::libs::pbit::cigar_delta::unpack_cigar;
use pgr::libs::pbit::collection::SegmentDesc;
use pgr::libs::pbit::decompressor::Decompressor;
use pgr::libs::pbit::format::DeltaEncoding;
use std::io::Write;

/// Build the clap subcommand for to-paf.
pub fn make_subcommand() -> Command {
    Command::new("to-paf")
        .about("Exports PAF alignments embedded in a pbit archive")
        .after_help(
            r###"
Exports the alignments embedded in a pbit archive as standard PAF (12
mandatory columns + `cg:Z` CIGAR).

Each PAF record corresponds to one CIGAR- or Identity-encoded segment: the
part of a sample contig that was aligned to the reference during
`pbit create --paf`. LZ-diff/Raw segments (uncovered parts) are not exported —
they carry no alignment information.

Notes:
* pbit files are binary and require random access (seeking)
* Does not support stdin or gzipped inputs

Examples:
1. Export all sample alignments:
   pgr pbit to-paf input.pbit -o out.paf

2. Export one sample's alignments:
   pgr pbit to-paf input.pbit -s sample1 -o out.paf
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input pbit file to process",
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(crate::cmd_pgr::args::pbit_sample_filter_arg(
            "Restrict output to this sample",
        ))
}

fn cigar_string(ops: &[CigarOp]) -> String {
    // Merge adjacent same-op runs (segment concatenation may split a run).
    let mut s = String::new();
    let mut prev_op: Option<char> = None;
    let mut prev_len: u32 = 0;
    for op in ops {
        let o = op.op();
        if Some(o) == prev_op {
            prev_len += op.len();
        } else {
            if let Some(p) = prev_op {
                s.push_str(&format!("{}{}", prev_len, p));
            }
            prev_op = Some(o);
            prev_len = op.len();
        }
    }
    if let Some(p) = prev_op {
        s.push_str(&format!("{}{}", prev_len, p));
    }
    s
}

fn cigar_matches(ops: &[CigarOp]) -> u32 {
    ops.iter()
        .map(|op| if op.op() == '=' { op.len() } else { 0 })
        .sum()
}

/// PAF block length: matches + mismatches + insertion + deletion bases.
fn cigar_block_length(ops: &[CigarOp]) -> u32 {
    ops.iter()
        .map(|op| {
            if matches!(op.op(), '=' | 'X' | 'I' | 'D') {
                op.len()
            } else {
                0
            }
        })
        .sum()
}

/// Rebuild a `cs:Z` string from a CIGAR, the X/I base stream, and the
/// reference interval, matching `maf to-paf`'s `cs_from_alignment` format:
/// `:run` for matches, `*RQ` for substitutions, `+Q` for insertions,
/// `-R` for deletions (uppercase).
fn cs_string(ops: &[CigarOp], xi: &[u8], ref_dna: &[u8]) -> String {
    let mut cs = String::new();
    let mut run = 0usize;
    let mut xi_pos = 0usize;
    let mut ref_pos = 0usize;
    let flush = |run: &mut usize, cs: &mut String| {
        if *run > 0 {
            cs.push(':');
            cs.push_str(&run.to_string());
            *run = 0;
        }
    };
    for op in ops {
        match op.op() {
            '=' => {
                run += op.len() as usize;
                ref_pos += op.len() as usize;
            }
            'X' => {
                flush(&mut run, &mut cs);
                for _ in 0..op.len() {
                    let q = xi[xi_pos];
                    xi_pos += 1;
                    let r = ref_dna[ref_pos];
                    ref_pos += 1;
                    cs.push('*');
                    cs.push(r.to_ascii_uppercase() as char);
                    cs.push(q.to_ascii_uppercase() as char);
                }
            }
            'I' => {
                flush(&mut run, &mut cs);
                for _ in 0..op.len() {
                    let q = xi[xi_pos];
                    xi_pos += 1;
                    cs.push('+');
                    cs.push(q.to_ascii_uppercase() as char);
                }
            }
            'D' => {
                flush(&mut run, &mut cs);
                for _ in 0..op.len() {
                    let r = ref_dna[ref_pos];
                    ref_pos += 1;
                    cs.push('-');
                    cs.push(r.to_ascii_uppercase() as char);
                }
            }
            _ => {}
        }
    }
    flush(&mut run, &mut cs);
    cs
}

/// Execute the to-paf command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infile = args
        .get_one::<String>("infile")
        .context("missing required argument: infile")?;
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let sample_filter = args.get_one::<String>("sample");

    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, vec![infile.as_str()])?;

    let mut dec = Decompressor::open(infile)
        .with_context(|| format!("Failed to open pbit file {}", infile))?;
    let mut writer = pgr::libs::io::writer(outfile)
        .with_context(|| format!("Failed to open writer for {}", outfile))?;

    let collection = dec.collection();
    let samples: Vec<String> = match sample_filter {
        Some(s) => {
            if !collection.samples.contains_key(s) {
                anyhow::bail!("sample '{}' not found in archive", s);
            }
            vec![s.clone()]
        }
        None => collection.samples.keys().cloned().collect(),
    };
    let mut work: Vec<(String, String, Vec<SegmentDesc>)> = Vec::new();
    for sample in &samples {
        if let Some(contigs) = collection.samples.get(sample) {
            for cs in contigs {
                work.push((sample.clone(), cs.contig_name.clone(), cs.segments.clone()));
            }
        }
    }
    let paf_data = dec.paf_data().to_vec();
    // ms lookup per sample: record_id → ms.
    let ms_of: std::collections::HashMap<String, std::collections::HashMap<u32, i32>> = paf_data
        .iter()
        .map(|(s, big, _)| (s.clone(), big.iter().cloned().collect()))
        .collect();
    let small_lines: std::collections::HashMap<String, Vec<String>> = paf_data
        .iter()
        .map(|(s, _, small)| (s.clone(), small.clone()))
        .collect();

    // Sample contig length table: (sample, contig) → total length.
    let mut contig_len: std::collections::HashMap<(String, String), u32> =
        std::collections::HashMap::new();
    for (sample, contig_name, segments) in &work {
        let mut total = 0u32;
        for seg in segments {
            let (meta, _) = dec.segment_payload(seg)?;
            total = total.max(seg.q_start + meta.raw_length);
        }
        contig_len.insert((sample.clone(), contig_name.clone()), total);
    }

    // Group CIGAR segments by source PAF record id (v1009) and rebuild each
    // big chain at chain level (segments merged by contiguity).
    let mut by_paf: std::collections::BTreeMap<u32, Vec<(String, String, SegmentDesc)>> =
        std::collections::BTreeMap::new();
    for (sample, contig_name, segments) in &work {
        for seg in segments {
            // Rebuild only complete big chains (present in the ms table);
            // incomplete big chains and small chains are output verbatim from
            // the stored PAF rows.
            if seg.paf_id != u32::MAX
                && ms_of
                    .get(sample)
                    .map(|m| m.contains_key(&seg.paf_id))
                    .unwrap_or(false)
            {
                by_paf.entry(seg.paf_id).or_default().push((
                    sample.clone(),
                    contig_name.clone(),
                    *seg,
                ));
            }
        }
    }
    for group in by_paf.values() {
        let mut group = group.clone();
        group.sort_by_key(|(_, _, seg)| seg.q_start);
        let (sample, contig_name, first) = group[0].clone();
        let (_, _, last) = group.last().cloned().unwrap();
        let (meta_first, _) = dec.segment_payload(&first)?;
        if !matches!(
            meta_first.encoding,
            DeltaEncoding::Cigar | DeltaEncoding::Identity
        ) {
            continue;
        }
        let mut ops = Vec::new();
        let mut xi_all: Vec<u8> = Vec::new();
        let mut matches = 0u32;
        let mut q_end = 0u32;
        // Concatenate segment CIGARs in chain direction: forward for '+',
        // reverse (RC) for '-'.
        let ordered: Vec<&(String, String, SegmentDesc)> = if meta_first.is_rev_comp {
            group.iter().rev().collect()
        } else {
            group.iter().collect()
        };
        for (_, _, seg) in ordered {
            let (meta, packed) = dec.segment_payload(seg)?;
            // Identity segments carry no payload; rebuild their CIGAR as a
            // single full-length '=' op so the chain merges losslessly.
            let (ops_seg, xi_seg) = match meta.encoding {
                DeltaEncoding::Cigar => {
                    unpack_cigar(&packed).with_context(|| "failed to unpack CIGAR".to_string())?
                }
                DeltaEncoding::Identity => {
                    (vec![CigarOp::try_new(meta.raw_length, '=')?], Vec::new())
                }
                _ => continue,
            };
            ops.extend_from_slice(&ops_seg);
            xi_all.extend_from_slice(&xi_seg);
            matches += cigar_matches(&ops_seg);
            q_end = q_end.max(seg.q_start + meta.raw_length);
        }
        let qstart = first.q_start;
        let strand = if meta_first.is_rev_comp { '-' } else { '+' };
        let Some((t_name, _, t_len)) = dec.ref_group_location(first.ref_group_id) else {
            continue;
        };
        let g_first = dec.ref_seg_starts()[first.ref_group_id as usize];
        let g_last = dec.ref_seg_starts()[last.ref_group_id as usize];
        let (loc_first, loc_last) = (
            dec.ref_group_location(first.ref_group_id),
            dec.ref_group_location(last.ref_group_id),
        );
        let (Some((_, t_seg_first, _)), Some((_, t_seg_last, _))) = (loc_first, loc_last) else {
            continue;
        };
        // '+' strand: CIGAR runs forward; segment ref_start maps to q_start.
        // '-' strand: CIGAR runs RC(query) vs forward(target); a segment's
        // ref_start maps to its q_end (target side) and ref_end to its
        // q_start. Chain tstart = target of the last segment's q_end
        // (= its ref_start); chain tend = target of the first segment's
        // q_start (= its ref_end).
        let (t_start, t_end) = if meta_first.is_rev_comp {
            (
                t_seg_last as u64 + (last.ref_start as u64 - g_last),
                t_seg_first as u64 + (first.ref_end as u64 - g_first),
            )
        } else {
            (
                t_seg_first as u64 + (first.ref_start as u64 - g_first),
                t_seg_last as u64 + (last.ref_end as u64 - g_last),
            )
        };
        let q_len = contig_len
            .get(&(sample.clone(), contig_name.clone()))
            .copied()
            .unwrap_or(q_end);
        let cg = cigar_string(&ops);
        let gi = gap_compressed_identity(&ops);
        let bi = block_identity(&ops);
        // Rebuild cs:Z from the concatenated CIGAR, X/I bases, and the
        // reference interval (global coordinates: '+' from the first
        // segment's ref_start to the last segment's ref_end; '-' the reverse).
        let (ref_g_start, ref_g_end) = if meta_first.is_rev_comp {
            (last.ref_start as u64, first.ref_end as u64)
        } else {
            (first.ref_start as u64, last.ref_end as u64)
        };
        let ref_dna = dec.read_ref_interval(ref_g_start, ref_g_end)?;
        let cs = cs_string(&ops, &xi_all, &ref_dna);
        let ms = ms_of
            .get(&sample)
            .and_then(|m| m.get(&first.paf_id))
            .copied();
        let mut line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255\tgi:f:{:.6}\tbi:f:{:.6}\tcg:Z:{}\tcs:Z:{}",
            contig_name,
            q_len,
            qstart,
            q_end,
            strand,
            t_name,
            t_len,
            t_start,
            t_end,
            matches,
            cigar_block_length(&ops),
            gi,
            bi,
            cg
            ,
            cs
        );
        if let Some(ms) = ms {
            line.push_str(&format!("\tms:i:{}", ms));
        }
        writeln!(writer, "{}", line)?;
    }
    // Verbatim small-chain rows.
    for sample in &samples {
        if let Some(lines) = small_lines.get(sample) {
            for line in lines {
                writeln!(writer, "{}", line)?;
            }
        }
    }
    Ok(())
}
