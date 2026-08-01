//! SD interval clustering and cluster FASTA extraction (BISER cluster, PAF-based).

use crate::libs::ds::Dsu;
use crate::libs::paf::parser::parse_paf;
use anyhow::Context;
use intspan::Range;
use std::collections::HashMap;
use std::io::{BufRead, Write};

/// One SD mate interval.
#[derive(Debug, Clone)]
pub struct SdInterval {
    /// FASTA header: `{species}#{chrom}{strand}#{start}#{end}` (0-based coords).
    pub name: String,
    /// Chromosome (species prefix stripped).
    pub chrom: String,
    pub strand: char,
    pub start: i32,
    pub end: i32,
    /// Sequence extracted from the genome (reverse-complemented if '-').
    pub seq: Vec<u8>,
}

/// A cluster of overlapping SD intervals.
pub struct SdCluster {
    pub intervals: Vec<SdInterval>,
}

/// Split a chainnet-style `{species}.{chrom}` name; falls back to `?` species.
fn split_species_name(full: &str) -> (String, String) {
    match full.split_once('.') {
        Some((sp, chr)) => (sp.to_string(), chr.to_string()),
        None => ("?".to_string(), full.to_string()),
    }
}

/// Cluster overlapping SD mates from a PAF, extract sequences from `genome`,
/// and write one FASTA per cluster into `outdir`. Returns the clusters.
///
/// Two mates of the same PAF record always share a cluster; intervals
/// overlapping on the same chromosome are unioned as well (BISER interval
/// coloring semantics, notes/references/biser.md §6.3.4).
pub fn cluster_paf<R: BufRead>(
    paf_reader: R,
    genome: &str,
    outdir: &str,
) -> anyhow::Result<Vec<SdCluster>> {
    let records = parse_paf(paf_reader)?;
    if records.is_empty() {
        return Ok(Vec::new());
    }

    let mut intervals: Vec<SdInterval> = Vec::new();
    let mut index: HashMap<(String, i32, i32), usize> = HashMap::new();
    let mut dsu = Dsu::new(records.len() * 2);

    let mut add = |chrom: &str, start: u32, end: u32, strand: char| -> usize {
        let key = (chrom.to_string(), start as i32, end as i32);
        if let Some(&i) = index.get(&key) {
            return i;
        }
        let (species, chr) = split_species_name(chrom);
        let i = intervals.len();
        intervals.push(SdInterval {
            name: format!("{species}#{chr}{strand}#{start}#{end}"),
            chrom: chr,
            strand,
            start: start as i32,
            end: end as i32,
            seq: Vec::new(),
        });
        index.insert(key, i);
        i
    };

    for rec in &records {
        let qi = add(&rec.query_name, rec.query_start, rec.query_end, rec.strand);
        let ti = add(&rec.target_name, rec.target_start, rec.target_end, '+');
        dsu.union(qi, ti); // both mates of one hit belong to the same cluster
    }

    // Union overlapping intervals on the same chromosome.
    let mut by_chrom: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, iv) in intervals.iter().enumerate() {
        by_chrom.entry(iv.chrom.clone()).or_default().push(i);
    }
    for idxs in by_chrom.values_mut() {
        idxs.sort_by_key(|&i| intervals[i].start);
        for w in idxs.windows(2) {
            let (a, b) = (w[0], w[1]);
            if intervals[b].start < intervals[a].end {
                dsu.union(a, b);
            }
        }
    }

    // Group intervals by connected component.
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..intervals.len() {
        groups.entry(dsu.find(i)).or_default().push(i);
    }

    // Extract sequences and write one FASTA per cluster.
    let (mut reader, loc_of) = crate::libs::loc::open_indexed(genome, false)?;
    let mut clusters: Vec<SdCluster> = Vec::new();
    for (_root, mut idxs) in groups {
        idxs.sort_by_key(|&i| (intervals[i].chrom.clone(), intervals[i].start));
        let mut members = Vec::new();
        for i in idxs {
            let iv = &mut intervals[i];
            let rg = Range::from_full(
                &iv.name,
                &iv.chrom,
                &iv.strand.to_string(),
                iv.start + 1,
                iv.end,
            );
            let seq = crate::libs::loc::fetch_range_seq(&mut reader, &loc_of, &rg)
                .with_context(|| format!("failed to fetch {}", iv.name))?;
            iv.seq = seq.into_bytes();
            members.push(iv.clone());
        }
        clusters.push(SdCluster { intervals: members });
    }
    clusters.sort_by_key(|c| c.intervals[0].chrom.clone());

    std::fs::create_dir_all(outdir)?;
    for (ci, cluster) in clusters.iter().enumerate() {
        let path = std::path::Path::new(outdir).join(format!("cluster_{}.fa", ci + 1));
        let mut w = std::io::BufWriter::new(std::fs::File::create(&path)?);
        for iv in &cluster.intervals {
            writeln!(w, ">{}", iv.name)?;
            for chunk in iv.seq.chunks(80) {
                writeln!(w, "{}", String::from_utf8_lossy(chunk))?;
            }
        }
    }
    Ok(clusters)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn species_name_split() {
        assert_eq!(
            split_species_name("mg1655.NC_000913"),
            ("mg1655".into(), "NC_000913".into())
        );
        assert_eq!(
            split_species_name("NC_000913"),
            ("?".into(), "NC_000913".into())
        );
    }
}
