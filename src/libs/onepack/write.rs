//! Write-side `.1aln` generation primitives.
//!
//! This is the P6 building block of the `.1aln` migration. It resamples a
//! base-level CIGAR into trace points ([`cigar2tp`], the only genuinely new
//! write-side algorithm) and writes the GDB skeleton plus alignment records
//! ([`write_skeleton_contigs`], [`write_aln_record`]) into a ONEcode container.
//! It mirrors FastGA's `PAFtoALN.c cigar2tp`, `GDB.c Write_Skeleton` and
//! `alncode.c Write_Aln_Overlap/Write_Aln_Trace`.
//!
//! The orchestration that maps a concrete source format (PAF-with-CIGAR, MAF)
//! onto these primitives lives in the CLI layer (`src/cmd_pgr`), not here.

use anyhow::Result;
use std::collections::HashMap;

use super::container::aln_schema_text;
use super::container::{Field, List, Writer};
use super::schema::parse_schema_text;
use crate::libs::paf::cigar::CigarOp;

/// The result of resampling a CIGAR into trace points.
pub struct TracePoints {
    /// `b`-coordinate deltas (`T` line), one per trace interval.
    pub tpoints: Vec<i64>,
    /// Per-interval difference counts (`X` line), same length as `tpoints`.
    pub tdiffs: Vec<i64>,
    /// Total differences (substitutions + indels) over the path.
    pub diffs: i64,
}

/// Resample a base-level CIGAR into trace points at `tspace` spacing.
///
/// The path starts at `(abpos, bbpos)` on the `a`/`b` axes. A trace point is
/// emitted each time the `a`-coordinate crosses a `tspace` multiple; each point
/// records the `b`-advance and the difference-advance since the previous point.
/// Mirrors FastGA `PAFtoALN.c cigar2tp`.
///
/// `=`/`M` advance both axes without accumulating differences, `X` advances
/// both and counts as differences, `I` advances only `a` (an insertion in `b`),
/// and `D` advances only `b` (a deletion in `b`). Deletions never cross an
/// `a`-boundary so they are folded into the next emitted point.
pub fn cigar2tp(ops: &[CigarOp], abpos: i64, bbpos: i64, tspace: i64) -> TracePoints {
    let mut apos = abpos;
    let mut bpos = bbpos;
    let mut diff: i64 = 0;
    let mut dlast: i64 = 0;
    let mut blast: i64 = bbpos;
    let mut anext = (abpos / tspace) * tspace + tspace;
    let mut tpoints: Vec<i64> = Vec::new();
    let mut tdiffs: Vec<i64> = Vec::new();

    // Emit a trace point at the current position.
    let emit = |tpoints: &mut Vec<i64>,
                tdiffs: &mut Vec<i64>,
                diff: i64,
                bpos: i64,
                dlast: &mut i64,
                blast: &mut i64| {
        tdiffs.push(diff - *dlast);
        tpoints.push(bpos - *blast);
        *blast = bpos;
        *dlast = diff;
    };

    for op in ops {
        let len = op.len() as i64;
        match op.op() {
            '=' | 'M' => {
                let mut remaining = len;
                while apos + remaining > anext {
                    let inc = anext - apos;
                    apos += inc;
                    bpos += inc;
                    remaining -= inc;
                    anext += tspace;
                    emit(
                        &mut tpoints,
                        &mut tdiffs,
                        diff,
                        bpos,
                        &mut dlast,
                        &mut blast,
                    );
                }
                apos += remaining;
                bpos += remaining;
            }
            'X' => {
                let mut remaining = len;
                while apos + remaining > anext {
                    let inc = anext - apos;
                    apos += inc;
                    bpos += inc;
                    diff += inc;
                    remaining -= inc;
                    anext += tspace;
                    emit(
                        &mut tpoints,
                        &mut tdiffs,
                        diff,
                        bpos,
                        &mut dlast,
                        &mut blast,
                    );
                }
                apos += remaining;
                bpos += remaining;
                diff += remaining;
            }
            'I' => {
                let mut remaining = len;
                while apos + remaining > anext {
                    let inc = anext - apos;
                    apos += inc;
                    diff += inc;
                    remaining -= inc;
                    anext += tspace;
                    emit(
                        &mut tpoints,
                        &mut tdiffs,
                        diff,
                        bpos,
                        &mut dlast,
                        &mut blast,
                    );
                }
                apos += remaining;
                diff += remaining;
            }
            // 'D': advance `b` and differences only; folded into the next point.
            'D' => {
                bpos += len;
                diff += len;
            }
            _ => {}
        }
    }
    // Emit the final partial interval.
    if apos > anext - tspace {
        emit(
            &mut tpoints,
            &mut tdiffs,
            diff,
            bpos,
            &mut dlast,
            &mut blast,
        );
    }
    TracePoints {
        tpoints,
        tdiffs,
        diffs: diff,
    }
}

/// Open a new `.1aln` writer using the embedded `aln` schema.
///
/// Writes the ASCII prolog (file type `aln`, provenance). The caller then adds
/// the source references (via [`Writer::add_reference`]) before writing any
/// line, emits the `t` line carrying `tspace` via [`write_tspace`], then the
/// skeleton(s) and records, before calling [`Writer::close`].
///
/// References must be added before the first [`Writer::write_line`], since the
/// ASCII header is flushed on the first write.
pub fn open_aln_writer(path: &str) -> Result<Writer> {
    let schema = parse_schema_text(aln_schema_text())?;
    let mut w = Writer::open(path, schema, true)?;
    w.add_provenance("pgr", env!("CARGO_PKG_VERSION"), "pgr 1aln to-1aln", "");
    Ok(w)
}

/// Write the `t` line carrying the trace-point spacing `tspace`.
pub fn write_tspace(w: &mut Writer, tspace: i64) -> Result<()> {
    w.write_line(b't', &[Field::Int(tspace)], None)
}

/// Write a GDB skeleton where each sequence is its own scaffold + contig.
///
/// `seqs` is `(name, length)` in the order the contigs should be indexed. This
/// is the contig-level case (no scaffold grouping): each name maps to a contig
/// index equal to its position in `seqs`. Returns that name → contig index map.
pub fn write_skeleton_contigs(
    w: &mut Writer,
    seqs: &[(String, i64)],
) -> Result<HashMap<String, usize>> {
    w.write_line(b'g', &[], None)?;
    let mut map = HashMap::new();
    for (i, (name, slen)) in seqs.iter().enumerate() {
        w.write_line(
            b'S',
            &[Field::Int(name.len() as i64)],
            Some(&List::Bytes(name.as_bytes().to_vec())),
        )?;
        w.write_line(b'C', &[Field::Int(*slen)], None)?;
        map.insert(name.clone(), i);
    }
    Ok(map)
}

/// Write one alignment record (`A` plus optional `R`, `D`, `T`, `X` lines).
///
/// `comp` marks the query (`b`) side as reverse-complemented. `tpoints`/`tdiffs`
/// are the decoded trace lists (see [`cigar2tp`]).
#[allow(clippy::too_many_arguments)] // mirrors the fixed 6-field `A` + `D`/`T`/`X` layout
pub fn write_aln_record(
    w: &mut Writer,
    aread: i64,
    abpos: i64,
    aepos: i64,
    bread: i64,
    bbpos: i64,
    bepos: i64,
    comp: bool,
    diffs: i64,
    tpoints: &[i64],
    tdiffs: &[i64],
) -> Result<()> {
    w.write_line(
        b'A',
        &[
            Field::Int(aread),
            Field::Int(abpos),
            Field::Int(aepos),
            Field::Int(bread),
            Field::Int(bbpos),
            Field::Int(bepos),
        ],
        None,
    )?;
    if comp {
        w.write_line(b'R', &[], None)?;
    }
    w.write_line(b'D', &[Field::Int(diffs)], None)?;
    w.write_line(
        b'T',
        &[Field::Int(tpoints.len() as i64)],
        Some(&List::Ints(tpoints.to_vec())),
    )?;
    w.write_line(
        b'X',
        &[Field::Int(tdiffs.len() as i64)],
        Some(&List::Ints(tdiffs.to_vec())),
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::container::Reader;
    use super::*;
    use crate::libs::paf::cigar::parse_cigar;

    fn temp_path(tag: &str) -> String {
        let mut p = std::env::temp_dir();
        p.push(format!("onecode_write_{tag}_{}.1aln", std::process::id()));
        p.to_str().unwrap().to_string()
    }

    #[test]
    fn cigar2tp_matches_fastga_trace_semantics() {
        // A short alignment (61 a-bases < tspace) yields a single trace point.
        // FastGA's point count is ((aepos-1)/tspace - abpos/tspace) + 1 = 1.
        let ops = parse_cigar("25=1X24=1I10=").unwrap();
        let tp = cigar2tp(&ops, 0, 0, 100);
        assert_eq!(tp.tpoints.len(), 1);
        // diffs = 1 X + 1 I.
        assert_eq!(tp.diffs, 2);
        // The single closing point records the full b-advance (60, the 1 I does
        // not advance b) and attributes both differences.
        assert_eq!(tp.tpoints[0], 60);
        assert_eq!(tp.tdiffs[0], 2);
    }

    #[test]
    fn cigar2tp_emits_point_on_boundary_crossing() {
        // A 120-bp a-span with abpos=0 crosses the 100 boundary once → two
        // trace points. The X sits exactly at apos=100, so the boundary point
        // records no b-advance/diff; the closing point records the rest.
        let ops = parse_cigar("100=1X19=").unwrap();
        let tp = cigar2tp(&ops, 0, 0, 100);
        assert_eq!(tp.tpoints.len(), 2);
        // The X sits at apos=100, so the boundary point records the 100 b-bases
        // already consumed and zero diffs; the closing point records the rest.
        assert_eq!(tp.tpoints[0], 100);
        assert_eq!(tp.tdiffs[0], 0);
        assert_eq!(tp.tpoints[1], 20);
        assert_eq!(tp.tdiffs[1], 1);
        assert_eq!(tp.diffs, 1);
    }

    #[test]
    fn write_read_round_trip() {
        let path = temp_path("roundtrip");
        let ops = parse_cigar("25=1X24=1I10=").unwrap();
        {
            let mut w = open_aln_writer(&path).unwrap();
            w.add_reference("ref.fa", 1);
            write_tspace(&mut w, 100).unwrap();
            let a_map = write_skeleton_contigs(&mut w, &[("r1".to_string(), 1000)]).unwrap();
            let b_map = write_skeleton_contigs(&mut w, &[("q1".to_string(), 1000)]).unwrap();
            assert_eq!(a_map["r1"], 0);
            assert_eq!(b_map["q1"], 0);
            let tp = cigar2tp(&ops, 0, 0, 100);
            write_aln_record(
                &mut w,
                0,
                0,
                60,
                0,
                0,
                59,
                false,
                tp.diffs,
                &tp.tpoints,
                &tp.tdiffs,
            )
            .unwrap();
            w.close().unwrap();
        }
        {
            let r = Reader::open(&path).unwrap();
            assert_eq!(r.references.len(), 1);
            assert_eq!(r.references[0].count, 1);
        }
        let _ = std::fs::remove_file(&path);
    }
}
