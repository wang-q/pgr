//! SD interval clustering and cluster FASTA extraction (BISER cluster, PAF-based).

use crate::libs::ds::Dsu;
use crate::libs::ds::Range;
use crate::libs::paf::parser::parse_paf;
use anyhow::Context;
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
/// coloring semantics, notes/design/sd.md §4.4).
pub fn cluster_paf<R: BufRead>(
    paf_reader: R,
    genome: &str,
    outdir: &str,
) -> anyhow::Result<Vec<SdCluster>> {
    let records = parse_paf(paf_reader)?;
    if records.is_empty() {
        return Ok(Vec::new());
    }
    // The `.loc` random-access index only supports plain and BGZF files; a
    // plain-gzip (non-BGZF) genome is decompressed into a temp file so the
    // cluster step accepts the same inputs as `sd search` / `sd align`.
    let (genome_path, _tmp) = plain_gz_to_temp(genome)?;

    let mut intervals: Vec<SdInterval> = Vec::new();
    let mut index: HashMap<(String, String, i32, i32, char), usize> = HashMap::new();
    let mut dsu = Dsu::new(records.len() * 2);

    let mut add = |chrom: &str, start: u32, end: u32, strand: char| -> usize {
        let (species, chr) = split_species_name(chrom);
        let key = (
            species.clone(),
            chr.clone(),
            start as i32,
            end as i32,
            strand,
        );
        if let Some(&i) = index.get(&key) {
            return i;
        }
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
        // Adjacent pairs alone miss nested intervals (A contains B and C,
        // B/C disjoint: A overlaps C but C never meets B in the scan).
        // Union each interval with the furthest-extending one seen so far:
        // any overlap with an earlier interval implies an overlap with it.
        let mut last = idxs[0];
        for &i in &idxs[1..] {
            if intervals[i].start < intervals[last].end {
                dsu.union(i, last);
            }
            if intervals[i].end > intervals[last].end {
                last = i;
            }
        }
    }

    // Group intervals by connected component, then order the clusters by
    // their first interval's (chrom, start). A HashMap iteration order is
    // randomized per process, so the cluster numbers (and `sd run`'s global
    // set_id renumbering) must not depend on it.
    let mut groups: Vec<Vec<usize>> = {
        let mut by_root: HashMap<usize, Vec<usize>> = HashMap::new();
        for i in 0..intervals.len() {
            by_root.entry(dsu.find(i)).or_default().push(i);
        }
        by_root.into_values().collect()
    };
    for idxs in &mut groups {
        idxs.sort_by_key(|&i| (intervals[i].chrom.clone(), intervals[i].start));
    }
    groups.sort_by_key(|idxs| {
        let first = idxs[0];
        (intervals[first].chrom.clone(), intervals[first].start)
    });

    // Extract sequences and write one FASTA per cluster.
    let (mut reader, loc_of) = crate::libs::loc::open_indexed(&genome_path, false)?;
    let mut clusters: Vec<SdCluster> = Vec::new();
    for idxs in groups {
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
    // Remove stale cluster files from a previous run into the same
    // directory; a leftover `cluster_2.fa` would otherwise be silently
    // consumed by `sd run` / `sd decompose` as if it were current output.
    // Only pgr's own `cluster_<N>.fa` naming is touched, never other files.
    if let Ok(entries) = std::fs::read_dir(outdir) {
        for entry in entries.flatten() {
            let name = entry.file_name().to_string_lossy().into_owned();
            let num = name
                .strip_prefix("cluster_")
                .and_then(|rest| rest.strip_suffix(".fa"));
            if num.is_some_and(|n| n.parse::<u32>().is_ok()) {
                let _ = std::fs::remove_file(entry.path());
            }
        }
    }
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

/// Return a path the `.loc` index can open: `genome` itself, or a temp
/// decompression when it is a plain-gzip (non-BGZF) `.gz` file.
fn plain_gz_to_temp(genome: &str) -> anyhow::Result<(String, Option<tempfile::TempDir>)> {
    let is_gz = std::path::Path::new(genome)
        .extension()
        .and_then(|e| e.to_str())
        == Some("gz");
    if !is_gz || crate::libs::io::is_bgzf(genome) {
        return Ok((genome.to_string(), None));
    }
    let dir = tempfile::TempDir::new()?;
    let out = dir.path().join("genome.plain.fa");
    let mut reader = crate::libs::io::reader(genome)?;
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&out)?);
    std::io::copy(&mut reader, &mut writer)?;
    Ok((out.to_string_lossy().into_owned(), Some(dir)))
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

    #[test]
    fn nested_overlapping_intervals_form_one_cluster() {
        // A [0,100] contains B [10,20] and C [30,40]; B and C do not overlap
        // each other, but both overlap A, so all three share a cluster.
        let dir = tempfile::tempdir().unwrap();
        let genome = dir.path().join("g.fa");
        std::fs::write(&genome, format!(">chr\n{}\n", "A".repeat(120))).unwrap();
        let paf = "\
chr\t120\t0\t100\t+\tchr\t120\t0\t100\t100\t100\t255\n\
chr\t120\t10\t20\t+\tchr\t120\t10\t20\t10\t10\t255\n\
chr\t120\t30\t40\t+\tchr\t120\t30\t40\t10\t10\t255\n";
        let outdir = dir.path().join("clusters");
        let clusters = cluster_paf(
            std::io::Cursor::new(paf.as_bytes()),
            genome.to_str().unwrap(),
            outdir.to_str().unwrap(),
        )
        .unwrap();
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].intervals.len(), 3);
    }

    #[test]
    fn same_coordinates_on_opposite_strands_are_distinct_copies() {
        // Regression: the interval dedup key used to ignore the strand (and
        // species), so two copies occupying the same coordinates on opposite
        // strands (an inverted self-repeat) were collapsed into one interval
        // with the first strand, losing the second copy.
        let dir = tempfile::tempdir().unwrap();
        let genome = dir.path().join("g.fa");
        std::fs::write(&genome, format!(">chr\n{}\n", "A".repeat(1000))).unwrap();
        // Record 1 aligns query [100,200] on the minus strand; record 2
        // aligns the same query interval on the plus strand.
        let paf = "\
chr\t1000\t100\t200\t-\tchr\t1000\t500\t600\t100\t100\t255\n\
chr\t1000\t100\t200\t+\tchr\t1000\t700\t800\t100\t100\t255\n";
        let outdir = dir.path().join("clusters");
        let clusters = cluster_paf(
            std::io::Cursor::new(paf.as_bytes()),
            genome.to_str().unwrap(),
            outdir.to_str().unwrap(),
        )
        .unwrap();
        assert_eq!(clusters.len(), 1);
        let same_coords: Vec<&SdInterval> = clusters[0]
            .intervals
            .iter()
            .filter(|iv| iv.chrom == "chr" && iv.start == 100 && iv.end == 200)
            .collect();
        assert_eq!(
            same_coords.len(),
            2,
            "plus and minus copies must both survive: {:?}",
            clusters[0].intervals
        );
        let strands: Vec<char> = same_coords.iter().map(|iv| iv.strand).collect();
        assert!(strands.contains(&'+') && strands.contains(&'-'));
    }

    #[test]
    fn stale_cluster_files_are_removed() {
        // Regression: re-running `sd cluster` into a directory that already
        // holds `cluster_N.fa` files from an earlier (larger) run left the
        // stale files behind, so `sd run` / `sd decompose` silently consumed
        // outdated families as if they were current output.
        let dir = tempfile::tempdir().unwrap();
        let genome = dir.path().join("g.fa");
        std::fs::write(&genome, format!(">chr\n{}\n", "A".repeat(300))).unwrap();
        let paf = "\
chr\t300\t0\t100\t+\tchr\t300\t100\t200\t100\t100\t255\n";
        let outdir = dir.path().join("clusters");
        std::fs::create_dir_all(&outdir).unwrap();
        std::fs::write(outdir.join("cluster_1.fa"), "stale").unwrap();
        std::fs::write(outdir.join("cluster_2.fa"), "stale").unwrap();
        std::fs::write(outdir.join("cluster_3.fa"), "stale").unwrap();
        // A file that is not pgr's cluster naming must survive.
        std::fs::write(outdir.join("notes.txt"), "keep me").unwrap();

        let clusters = cluster_paf(
            std::io::Cursor::new(paf.as_bytes()),
            genome.to_str().unwrap(),
            outdir.to_str().unwrap(),
        )
        .unwrap();
        assert_eq!(clusters.len(), 1);
        let mut names: Vec<String> = std::fs::read_dir(&outdir)
            .unwrap()
            .filter_map(|e| e.ok())
            .map(|e| e.file_name().to_string_lossy().into_owned())
            .collect();
        names.sort();
        assert_eq!(
            names,
            vec!["cluster_1.fa".to_string(), "notes.txt".to_string()],
            "stale cluster files must be removed, unrelated files kept"
        );
    }
}
