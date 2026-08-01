//! Core duplicon identification: greedy set cover over elementary SDs.

use crate::libs::paf::parser::parse_paf;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Write};

/// One elementary SD row (genome coordinates, from `pgr sd decompose`).
#[derive(Debug, Clone)]
pub struct ElemRow {
    pub species: String,
    pub chrom: String,
    pub begin: usize,
    pub end: usize,
    pub set_id: u32,
    pub length: usize,
    pub score: u32,
    pub strand: char,
}

/// Parse elementary SD BED rows (8 columns, decompose output).
pub fn read_elems<R: BufRead>(reader: R) -> anyhow::Result<Vec<ElemRow>> {
    let mut rows = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            continue;
        }
        rows.push(ElemRow {
            species: f[0].to_string(),
            chrom: f[1].to_string(),
            begin: f[2].parse()?,
            end: f[3].parse()?,
            set_id: f[4].parse()?,
            length: f[5].parse()?,
            score: f[6].parse()?,
            strand: f[7].chars().next().unwrap_or('+'),
        });
    }
    Ok(rows)
}

fn overlaps(a_begin: usize, a_end: usize, b_begin: usize, b_end: usize) -> bool {
    a_begin < b_end && b_begin < a_end
}

/// Strip a chainnet `{species}.{chrom}` prefix for chromosome comparison.
fn chrom_of(name: &str) -> &str {
    name.split_once('.').map(|(_, c)| c).unwrap_or(name)
}

/// Mark elementary SD sets that greedily cover all SD hits.
///
/// An elementary SD set (same `set_id`) covers an SD hit if any of its copies
/// overlaps the hit's query or target interval on the same chromosome.
/// Returns `true` per input row for rows whose set is selected as CORE.
pub fn cover_hits(rows: &[ElemRow], hits: &[crate::libs::paf::record::PafRecord]) -> Vec<bool> {
    if rows.is_empty() {
        return Vec::new();
    }

    // Group elementary copies by set_id.
    let mut set_order: Vec<u32> = Vec::new();
    let mut by_set: HashMap<u32, Vec<usize>> = HashMap::new();
    for (i, r) in rows.iter().enumerate() {
        if !by_set.contains_key(&r.set_id) {
            set_order.push(r.set_id);
        }
        by_set.entry(r.set_id).or_default().push(i);
    }

    // Elementary set -> covered SD hit indices.
    let mut coverage: HashMap<u32, HashSet<usize>> = HashMap::new();
    for (si, hit) in hits.iter().enumerate() {
        for &sid in &set_order {
            let covered = by_set[&sid].iter().any(|&ri| {
                let r = &rows[ri];
                (r.chrom == chrom_of(&hit.query_name) || r.chrom == chrom_of(&hit.target_name))
                    && ((r.chrom == chrom_of(&hit.query_name)
                        && overlaps(
                            r.begin,
                            r.end,
                            hit.query_start as usize,
                            hit.query_end as usize,
                        ))
                        || (r.chrom == chrom_of(&hit.target_name)
                            && overlaps(
                                r.begin,
                                r.end,
                                hit.target_start as usize,
                                hit.target_end as usize,
                            )))
            });
            if covered {
                coverage.entry(sid).or_default().insert(si);
            }
        }
    }

    // Greedy set cover over uncovered hits.
    let mut uncovered: HashSet<usize> = (0..hits.len()).collect();
    let mut core: HashSet<u32> = HashSet::new();
    while !uncovered.is_empty() {
        let best = set_order
            .iter()
            .filter(|sid| !core.contains(*sid))
            .map(|sid| {
                let c = coverage
                    .get(sid)
                    .map(|s| s.intersection(&uncovered).count())
                    .unwrap_or(0);
                (*sid, c)
            })
            .max_by_key(|&(_, c)| c);
        match best {
            Some((_, 0)) => break, // nothing covers the rest; stop
            Some((sid, _)) => {
                core.insert(sid);
                if let Some(c) = coverage.get(&sid) {
                    uncovered.retain(|u| !c.contains(u));
                }
            }
            None => break,
        }
    }

    rows.iter().map(|r| core.contains(&r.set_id)).collect()
}

/// Write elementary rows with a `CORE`/`non-core` marker column.
pub fn write_covered<W: Write>(
    rows: &[ElemRow],
    is_core: &[bool],
    writer: &mut W,
) -> anyhow::Result<()> {
    for (r, core) in rows.iter().zip(is_core) {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.species,
            r.chrom,
            r.begin,
            r.end,
            r.set_id,
            r.length,
            r.score,
            r.strand,
            if *core { "CORE" } else { "non-core" }
        )?;
    }
    Ok(())
}

/// Load SD hits and elementary rows, mark CORE, and write the annotated BED.
pub fn run_cover<R1: BufRead, R2: BufRead, W: Write>(
    hits_reader: R1,
    elems_reader: R2,
    writer: &mut W,
) -> anyhow::Result<()> {
    let hits = parse_paf(hits_reader)?;
    let rows = read_elems(elems_reader)?;
    let is_core = cover_hits(&rows, &hits);
    write_covered(&rows, &is_core, writer)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::paf::record::PafRecord;

    fn row(chrom: &str, b: usize, e: usize, sid: u32) -> ElemRow {
        ElemRow {
            species: "s".into(),
            chrom: chrom.into(),
            begin: b,
            end: e,
            set_id: sid,
            length: e - b,
            score: 1,
            strand: '+',
        }
    }

    fn hit(qb: u32, qe: u32, tb: u32, te: u32) -> PafRecord {
        PafRecord {
            query_name: "chr".into(),
            query_length: 1000,
            query_start: qb,
            query_end: qe,
            strand: '+',
            target_name: "chr".into(),
            target_length: 1000,
            target_start: tb,
            target_end: te,
            matches: 10,
            block_length: 10,
            mapq: 255,
            tags: vec![],
        }
    }

    #[test]
    fn greedy_cover_selects_minimal_sets() {
        // Elementary A covers hits 1-2, B covers hit 3 only; A is core.
        let rows = vec![row("chr", 0, 100, 1), row("chr", 500, 600, 2)];
        let hits = vec![hit(10, 50, 0, 0), hit(20, 60, 0, 0), hit(520, 550, 0, 0)];
        let core = cover_hits(&rows, &hits);
        assert_eq!(core, vec![true, true]);
    }

    #[test]
    fn no_overlap_no_cover() {
        let rows = vec![row("chr", 0, 100, 1)];
        let hits = vec![hit(500, 550, 0, 0)];
        let core = cover_hits(&rows, &hits);
        assert_eq!(core, vec![false]);
    }
}
