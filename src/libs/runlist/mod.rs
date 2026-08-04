//! Runlist interval operations, migrated from the external `spanr` CLI of the
//! intspan project.
//!
//! Operates on two input forms:
//! * `.rg` lines (`chr:start-end`, 1-based inclusive) for `cover`/`coverage`;
//! * runlist JSON (`{"chr": "start-end,..."}` or the multi form
//!   `{"name": {"chr": "..."}}`) for `span`/`compare`/`merge`.
//!
//! `coverage` computes per-position depth with a sweep line over sorted
//! start/end events (O(n log n)), which is the standard efficient approach
//! for depth aggregation; interval trees are not needed for pure depth.

use crate::libs::ds::IntSpan;
use anyhow::Context;
use coitrees::{BasicCOITree, Interval as CoiInterval, IntervalTree};
use std::collections::BTreeMap;
use std::io::BufRead;

/// Whether `range` is a usable `.rg` range: valid coordinates in ascending
/// order and within the representable maximum (`POS_INF - 1`).
pub fn usable_range(range: &crate::libs::ds::Range) -> bool {
    let max_coord = IntSpan::new().get_pos_inf();
    range.is_valid()
        && range.start() <= range.end()
        && range.start() <= &max_coord
        && range.end() <= &max_coord
}

/// Parse `.rg` lines (`chr:start-end`, 1-based inclusive) from `reader` into
/// a per-chromosome merged `IntSpan`. Lines starting with `#` and lines that
/// do not parse as valid ranges are skipped.
pub fn rg_to_set<R: BufRead>(reader: R) -> anyhow::Result<BTreeMap<String, IntSpan>> {
    let mut set: BTreeMap<String, IntSpan> = BTreeMap::new();
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let range = crate::libs::ds::Range::from_str(line);
        if !usable_range(&range) {
            continue;
        }
        set.entry(range.chr().clone())
            .or_default()
            .add_pair(*range.start(), *range.end());
    }
    Ok(set)
}

/// Parse `.rg` lines into per-chromosome half-open `[start, end+1)` interval
/// lists, preserving multiplicity for depth computation. Lines starting with
/// `#` and unparseable lines are skipped.
pub fn rg_to_intervals<R: BufRead>(reader: R) -> anyhow::Result<BTreeMap<String, Vec<(u32, u32)>>> {
    let mut iv_of: BTreeMap<String, Vec<(u32, u32)>> = BTreeMap::new();
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let range = crate::libs::ds::Range::from_str(line);
        if !usable_range(&range) {
            continue;
        }
        iv_of
            .entry(range.chr().clone())
            .or_default()
            .push((*range.start() as u32, *range.end() as u32 + 1));
    }
    Ok(iv_of)
}

/// Per-chromosome COITree index over `.rg` intervals, for overlap counting.
pub struct RgIndex {
    trees: BTreeMap<String, BasicCOITree<bool, u32>>,
}

impl RgIndex {
    /// Build the index from one or more `.rg` files; lines that do not parse
    /// as valid ranges are skipped.
    pub fn from_files(files: &[String]) -> anyhow::Result<Self> {
        let mut intervals_of: BTreeMap<String, Vec<CoiInterval<bool>>> = BTreeMap::new();
        for f in files {
            let reader = crate::reader(f)?;
            for line in reader.lines() {
                let line = line?;
                let range = crate::libs::ds::Range::from_str(&line);
                if !usable_range(&range) {
                    continue;
                }
                intervals_of
                    .entry(range.chr().clone())
                    .or_default()
                    .push(CoiInterval::new(*range.start(), *range.end(), true));
            }
        }
        let trees = intervals_of
            .into_iter()
            .map(|(chr, ivs)| (chr, BasicCOITree::new(&ivs)))
            .collect();
        Ok(Self { trees })
    }

    /// Number of intervals overlapping the inclusive range `[start, end]` on
    /// `chr` (0 when the chromosome has no indexed intervals).
    pub fn count(&self, chr: &str, start: i32, end: i32) -> usize {
        match self.trees.get(chr) {
            Some(tree) => {
                let mut n = 0usize;
                tree.query(start, end, |_| n += 1);
                n
            }
            None => 0,
        }
    }
}

/// Sorted, disjoint per-chromosome span lists for O(log n + k) intersection
/// queries against a runlist set.
pub struct SpanIndex {
    spans_of: BTreeMap<String, Vec<(i32, i32)>>,
}

impl SpanIndex {
    /// Build the index from a runlist set; each chromosome's spans are
    /// disjoint and sorted in ascending order (an `IntSpan` invariant).
    pub fn from_set(set: &BTreeMap<String, IntSpan>) -> Self {
        let spans_of = set
            .iter()
            .map(|(chr, ints)| (chr.clone(), ints.spans()))
            .collect();
        Self { spans_of }
    }

    /// Size of the intersection of `[start, end]` with the indexed runlist
    /// and the length of `[start, end]` (0 and `end - start + 1` when the
    /// chromosome is absent).
    pub fn overlap(&self, chr: &str, start: i32, end: i32) -> (i32, i32) {
        let length = end - start + 1;
        let size = match self.spans_of.get(chr) {
            Some(spans) => {
                // Spans are disjoint and sorted, so the ones overlapping
                // `[start, end]` form a contiguous range: those with
                // `end >= start` and `start <= end`.
                let first = spans.partition_point(|&(_, sp_end)| sp_end < start);
                let last = spans.partition_point(|&(sp_start, _)| sp_start <= end);
                let mut total: i64 = 0;
                for &(sp_start, sp_end) in &spans[first..last] {
                    total += i64::from(end.min(sp_end) - start.max(sp_start) + 1);
                }
                total.min(i64::from(i32::MAX)) as i32
            }
            None => 0,
        };
        (size, length)
    }
}

/// Proportion of `[start, end]` covered by `index`, plus the intersection
/// size and the range length (0/0.0 when the chromosome is absent).
pub fn range_prop(index: &SpanIndex, chr: &str, start: i32, end: i32) -> (f32, i32, i32) {
    let (size, length) = index.overlap(chr, start, end);
    let prop = if length == 0 {
        0.0
    } else {
        size as f32 / length as f32
    };
    (prop, length, size)
}

/// Regions of half-open `[start, end)` intervals covered by at least
/// `min_depth` intervals (sweep-line depth). Returns an inclusive `IntSpan`
/// runlist. `min_depth` of 0 is treated as 1 (only covered positions exist).
pub fn depth_at_least(ivs: &[(u32, u32)], min_depth: u32) -> IntSpan {
    depth_runs(ivs, min_depth).0
}

/// Per-depth regions (key = depth string) covered by at least `min_depth`
/// intervals (detailed mode of the spanr `coverage` command).
pub fn depth_by_level(ivs: &[(u32, u32)], min_depth: u32) -> BTreeMap<String, IntSpan> {
    depth_runs(ivs, min_depth).1
}

/// Sweep over sorted start/end events once, emitting either the merged
/// `>= min_depth` runlist plus runs grouped by their exact depth.
fn depth_runs(ivs: &[(u32, u32)], min_depth: u32) -> (IntSpan, BTreeMap<String, IntSpan>) {
    let min_depth = min_depth.max(1) as i64;
    let mut events: Vec<(i64, i64)> = Vec::with_capacity(ivs.len() * 2);
    for &(s, e) in ivs {
        events.push((s as i64, 1));
        events.push((e as i64, -1));
    }
    events.sort_unstable();

    let mut by_level: BTreeMap<String, IntSpan> = BTreeMap::new();
    let mut at_least = IntSpan::new();
    let mut depth = 0i64;
    let mut run_start: Option<i64> = None;
    let mut run_depth = 0i64;
    let mut i = 0usize;
    while i < events.len() {
        let pos = events[i].0;
        let mut delta = 0i64;
        while i < events.len() && events[i].0 == pos {
            delta += events[i].1;
            i += 1;
        }
        // Close the previous run (started at `run_start` with `run_depth`).
        if let Some(s) = run_start.take() {
            if pos > s {
                if run_depth >= min_depth {
                    by_level
                        .entry(run_depth.to_string())
                        .or_default()
                        .add_pair(s as i32, pos.saturating_sub(1) as i32);
                }
                if run_depth >= min_depth {
                    at_least.add_pair(s as i32, pos.saturating_sub(1) as i32);
                }
            }
        }
        depth += delta;
        run_start = Some(pos);
        run_depth = depth;
    }
    // A trailing run can only remain when the last event is a start (an
    // open-ended interval); close it at the last event position.
    if let Some(s) = run_start {
        let e = events.last().map(|x| x.0).unwrap_or(s);
        if e > s {
            if run_depth >= min_depth {
                by_level
                    .entry(run_depth.to_string())
                    .or_default()
                    .add_pair(s as i32, e.saturating_sub(1) as i32);
            }
            if run_depth >= min_depth {
                at_least.add_pair(s as i32, e.saturating_sub(1) as i32);
            }
        }
    }
    (at_least, by_level)
}

/// Operations of the `span` subcommand, applied per chromosome.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpanOp {
    Cover,
    Holes,
    Trim,
    Pad,
    Excise,
    Fill,
}

/// Apply `op` (with parameter `n` for trim/pad/excise/fill) to every
/// chromosome of `set`.
pub fn span_op(set: &BTreeMap<String, IntSpan>, op: SpanOp, n: i32) -> BTreeMap<String, IntSpan> {
    set.iter()
        .map(|(chr, is)| {
            let out = match op {
                SpanOp::Cover => is.cover(),
                SpanOp::Holes => is.holes(),
                SpanOp::Trim => is.trim(n),
                SpanOp::Pad => is.pad(n),
                SpanOp::Excise => is.excise(n),
                SpanOp::Fill => is.fill(n),
            };
            (chr.clone(), out)
        })
        .collect()
}

/// Set-comparison operations of the `compare` subcommand.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompareOp {
    Intersect,
    Union,
    Diff,
    Xor,
}

/// Compare the (possibly multi) `first` set against one or more single
/// `others`, filling missing chromosomes with empty sets.
pub fn compare_sets(
    first: &BTreeMap<String, BTreeMap<String, IntSpan>>,
    others: &[BTreeMap<String, IntSpan>],
    op: CompareOp,
) -> BTreeMap<String, BTreeMap<String, IntSpan>> {
    let mut chrs: Vec<String> = first
        .values()
        .flat_map(|s| s.keys().cloned())
        .chain(others.iter().flat_map(|s| s.keys().cloned()))
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    chrs.sort_unstable();

    let mut res: BTreeMap<String, BTreeMap<String, IntSpan>> = BTreeMap::new();
    for (name, s1) in first {
        let mut out: BTreeMap<String, IntSpan> = BTreeMap::new();
        for chr in &chrs {
            let mut acc = s1.get(chr).cloned().unwrap_or_default();
            for s2 in others {
                let other = s2.get(chr).cloned().unwrap_or_default();
                acc = match op {
                    CompareOp::Intersect => acc.intersect(&other),
                    CompareOp::Union => acc.union(&other),
                    CompareOp::Diff => acc.diff(&other),
                    CompareOp::Xor => acc.xor(&other),
                };
            }
            out.insert(chr.clone(), acc);
        }
        res.insert(name.clone(), out);
    }
    res
}

/// Read a runlist JSON file into a value map.
pub fn read_json(path: &str) -> anyhow::Result<BTreeMap<String, serde_json::Value>> {
    let mut reader = crate::reader(path)?;
    let mut s = String::new();
    std::io::Read::read_to_string(&mut reader, &mut s)
        .with_context(|| format!("failed to read runlist JSON {}", path))?;
    serde_json::from_str(&s).with_context(|| format!("failed to parse runlist JSON {}", path))
}

/// Convert a single runlist JSON map into per-chromosome `IntSpan`s.
///
/// Non-string values and unparseable runlist strings are errors rather than
/// being silently dropped (a multi runlist passed where a single one is
/// expected used to become an empty set).
pub fn json_to_set(
    json: &BTreeMap<String, serde_json::Value>,
) -> anyhow::Result<BTreeMap<String, IntSpan>> {
    let mut set = BTreeMap::new();
    for (chr, v) in json {
        let s = v
            .as_str()
            .ok_or_else(|| anyhow::anyhow!("runlist value for {} is not a string", chr))?;
        set.insert(
            chr.clone(),
            IntSpan::try_from(s).with_context(|| format!("invalid runlist for {}", chr))?,
        );
    }
    Ok(set)
}

/// Convert a runlist JSON map into `name -> chromosome -> IntSpan`; detects
/// the multi form (`{"name": {"chr": "..."}}`) automatically. An empty or
/// flat input is treated as a single set under the `__single__` key.
pub fn json_to_sets(
    json: &BTreeMap<String, serde_json::Value>,
) -> anyhow::Result<BTreeMap<String, BTreeMap<String, IntSpan>>> {
    let is_multi = json.values().next().map(|v| v.is_object()).unwrap_or(false);
    if is_multi {
        json.iter()
            .map(|(name, v)| {
                let inner = v.as_object().ok_or_else(|| {
                    anyhow::anyhow!("runlist value for {} is not an object", name)
                })?;
                let set =
                    json_to_set(&inner.iter().map(|(k, v)| (k.clone(), v.clone())).collect())?;
                Ok((name.clone(), set))
            })
            .collect()
    } else {
        let mut m = BTreeMap::new();
        m.insert("__single__".to_string(), json_to_set(json)?);
        Ok(m)
    }
}

/// Write the `name -> chromosome -> IntSpan` map as multi runlist JSON, or
/// as flat runlist JSON when it holds only the `__single__` set.
pub fn write_sets(
    output: &str,
    set_of: &BTreeMap<String, BTreeMap<String, IntSpan>>,
) -> anyhow::Result<()> {
    if set_of.len() == 1 && set_of.contains_key("__single__") {
        let json = crate::libs::ds::intspan::set2json(&set_of["__single__"]);
        crate::libs::ds::intspan::write_json(output, &json)
    } else {
        let json = crate::libs::ds::intspan::set2json_m(set_of);
        crate::libs::ds::intspan::write_json(output, &json)
    }
}

/// Merge several runlist JSON files into a multi runlist keyed by file stem
/// (`all = false` keeps only the first dot-separated segment).
pub fn merge_files(
    files: &[String],
    all: bool,
) -> anyhow::Result<BTreeMap<String, serde_json::Value>> {
    let mut out: BTreeMap<String, serde_json::Value> = BTreeMap::new();
    for f in files {
        let stem = std::path::Path::new(f)
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| anyhow::anyhow!("invalid file stem for {}", f))?;
        let key = if all {
            stem.to_string()
        } else {
            stem.split('.').next().unwrap_or(stem).to_string()
        };
        let json = read_json(f)?;
        out.insert(key, serde_json::to_value(json)?);
    }
    Ok(out)
}

/// Combine all sets of a multi runlist into one, applying `op` between the
/// first set and each subsequent one (spanr `combine` semantics).
pub fn combine_sets(
    set_of: &BTreeMap<String, BTreeMap<String, IntSpan>>,
    op: CompareOp,
) -> BTreeMap<String, IntSpan> {
    let mut names: Vec<String> = set_of.keys().cloned().collect();
    if names.is_empty() {
        return BTreeMap::new();
    }
    let first = names.remove(0);
    let mut wrapped = BTreeMap::new();
    wrapped.insert(first.clone(), set_of[&first].clone());
    let others: Vec<BTreeMap<String, IntSpan>> = names.iter().map(|n| set_of[n].clone()).collect();
    compare_sets(&wrapped, &others, op)
        .into_values()
        .next()
        .unwrap_or_default()
}

/// Ranges of one runlist set as `chr:start-end` lines (one per sub-span), or
/// only the longest span per chromosome with `longest`.
pub fn convert_set(set: &BTreeMap<String, IntSpan>, longest: bool) -> Vec<String> {
    let mut out = Vec::new();
    for (chr, ints) in set {
        let mut intses = ints.intses();
        if longest {
            if intses.is_empty() {
                continue;
            }
            intses.sort_by_cached_key(|e| -e.size());
            out.push(format!("{}:{}", chr, intses.first().unwrap()));
        } else {
            for sub in &intses {
                out.push(format!("{}:{}", chr, sub));
            }
        }
    }
    out
}

/// A runlist set covering every chromosome of `sizes` in full (1..size).
///
/// Non-positive sizes are rejected instead of panicking in
/// `IntSpan::from_pair` (a sizes file can contain 0-length contigs).
pub fn genome_set(sizes: &BTreeMap<String, i32>) -> anyhow::Result<BTreeMap<String, IntSpan>> {
    let mut set = BTreeMap::new();
    for (k, &v) in sizes {
        if v <= 0 {
            anyhow::bail!("invalid chromosome size {} for {}", v, k);
        }
        if v > IntSpan::new().get_pos_inf() {
            anyhow::bail!("chromosome size {} out of range for {}", v, k);
        }
        set.insert(k.clone(), IntSpan::from_pair(1, v));
    }
    Ok(set)
}

/// Subset of a runlist JSON whose top-level keys are in `names`.
pub fn some_json(
    json: &BTreeMap<String, serde_json::Value>,
    names: &std::collections::BTreeSet<String>,
) -> BTreeMap<String, serde_json::Value> {
    json.iter()
        .filter(|(k, _)| names.contains(*k))
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect()
}

/// Split a multi runlist JSON into `(key, compact JSON string)` pairs; errors
/// when a value is not itself an object.
pub fn split_json(
    json: &BTreeMap<String, serde_json::Value>,
) -> anyhow::Result<Vec<(String, String)>> {
    json.iter()
        .map(|(k, v)| {
            if !v.is_object() {
                anyhow::bail!("not a valid multi-key runlist json file (key {k})");
            }
            Ok((k.clone(), serde_json::to_string(v)?))
        })
        .collect()
}

/// Per-chromosome coverage stats as CSV lines (spanr `stat`); errors when a
/// chromosome of `set` is missing from `sizes`.
pub fn stat_lines(
    set: &BTreeMap<String, IntSpan>,
    sizes: &BTreeMap<String, i32>,
    all: bool,
    prefix: Option<&str>,
) -> anyhow::Result<String> {
    let mut lines = String::new();
    let mut all_length: i64 = 0;
    let mut all_size: i64 = 0;
    for chr in set.keys() {
        let length = *sizes
            .get(chr)
            .ok_or_else(|| anyhow::anyhow!("chromosome {chr} not found in sizes"))?;
        let size = set[chr].cardinality();
        if let Some(s) = prefix {
            lines.push_str(&format!("{},", s));
        }
        lines.push_str(&format!(
            "{},{},{},{:.4}\n",
            chr,
            length,
            size,
            size as f32 / length as f32
        ));
        all_length += i64::from(length);
        all_size += i64::from(size);
    }
    let mut all_line = format!(
        "{},{},{},{:.4}\n",
        "all",
        all_length,
        all_size,
        all_size as f64 / all_length as f64
    );
    if all {
        lines = String::new();
        all_line = all_line.replacen("all,", "", 1);
    }
    if let Some(s) = prefix {
        all_line.insert_str(0, &format!("{},", s));
    }
    lines.push_str(all_line.trim_end());
    Ok(lines)
}

/// Cross-set coverage stats as CSV lines (spanr `statop`); `set_op` is the
/// per-chromosome result of `op` between `s1` and `s2`.
#[allow(clippy::too_many_arguments)]
pub fn statop_lines(
    s1: &BTreeMap<String, IntSpan>,
    sizes: &BTreeMap<String, i32>,
    s2: &BTreeMap<String, IntSpan>,
    set_op: &BTreeMap<String, IntSpan>,
    all: bool,
    prefix: Option<&str>,
) -> anyhow::Result<String> {
    let mut lines = String::new();
    let mut all_length: i64 = 0;
    let mut all_size: i64 = 0;
    let mut all_s2_length: i64 = 0;
    let mut all_s2_size: i64 = 0;
    for chr in s1.keys() {
        let length = *sizes
            .get(chr)
            .ok_or_else(|| anyhow::anyhow!("chromosome {chr} not found in sizes"))?;
        let size = s1[chr].cardinality();
        // Missing chromosomes in `s2`/`set_op` count as empty (0).
        let s2_length = s2.get(chr).map(IntSpan::cardinality).unwrap_or(0);
        let s2_size = set_op.get(chr).map(IntSpan::cardinality).unwrap_or(0);
        let c1 = size as f64 / length as f64;
        let c2 = if s2_length == 0 {
            0.0
        } else {
            s2_size as f64 / s2_length as f64
        };
        let ratio = if (c1 - 0.0).abs() < 0.00001 {
            0.0
        } else {
            c2 / c1
        };
        if let Some(s) = prefix {
            lines.push_str(&format!("{},", s));
        }
        lines.push_str(&format!(
            "{},{},{},{},{},{:.4},{:.4},{:.4}\n",
            chr, length, size, s2_length, s2_size, c1, c2, ratio
        ));
        all_length += i64::from(length);
        all_size += i64::from(size);
        all_s2_length += i64::from(s2_length);
        all_s2_size += i64::from(s2_size);
    }
    let all_c1 = all_size as f64 / all_length as f64;
    let all_c2 = if all_s2_length == 0 {
        0.0
    } else {
        all_s2_size as f64 / all_s2_length as f64
    };
    let all_ratio = if (all_c1 - 0.0).abs() < 0.00001 {
        0.0
    } else {
        all_c2 / all_c1
    };
    let mut all_line = format!(
        "{},{},{},{},{},{:.4},{:.4},{:.4}\n",
        "all", all_length, all_size, all_s2_length, all_s2_size, all_c1, all_c2, all_ratio
    );
    if all {
        lines = String::new();
        all_line = all_line.replacen("all,", "", 1);
    }
    if let Some(s) = prefix {
        all_line.insert_str(0, &format!("{},", s));
    }
    lines.push_str(all_line.trim_end());
    Ok(lines)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn set(runs: &str) -> IntSpan {
        IntSpan::from(runs)
    }

    #[test]
    fn rg_to_set_merges_ranges() {
        let rg = "chr1:1-10\nchr1:5-15\nchr2(+):100-200\nbad line\n#comment\n";
        let s = rg_to_set(std::io::Cursor::new(rg)).unwrap();
        assert_eq!(s["chr1"].to_string(), "1-15");
        assert_eq!(s["chr2"].to_string(), "100-200");
        assert_eq!(s.len(), 2);
    }

    #[test]
    fn rg_to_set_skips_reversed_ranges() {
        // `chr1:10-5` parses as start > end and used to panic in add_pair.
        let s = rg_to_set(std::io::Cursor::new("chr1:10-5\nchr1:1-10\n")).unwrap();
        assert_eq!(s["chr1"].to_string(), "1-10");
        let ivs = rg_to_intervals(std::io::Cursor::new("chr1:10-5\nchr1:1-10\n")).unwrap();
        assert_eq!(ivs["chr1"], vec![(1, 11)]);

        // Coordinates above POS_INF - 1 are skipped instead of overflowing
        // `add_pair` (which stores upper + 1 as an edge).
        let s = rg_to_set(std::io::Cursor::new(
            "chr1:2147483647-2147483647\nchr1:1-10\n",
        ))
        .unwrap();
        assert_eq!(s["chr1"].to_string(), "1-10");
        let ivs = rg_to_intervals(std::io::Cursor::new("chr1:2147483647-2147483647\n")).unwrap();
        assert!(ivs.is_empty());
    }

    #[test]
    fn span_index_overlap() {
        let mut runset = BTreeMap::new();
        runset.insert("chr1".to_string(), set("1-10,20-30"));
        let idx = SpanIndex::from_set(&runset);
        // Partial overlaps at both ends of the query.
        assert_eq!(idx.overlap("chr1", 5, 25), (12, 21));
        // Query entirely in a gap.
        assert_eq!(idx.overlap("chr1", 11, 19), (0, 9));
        // Query touching span boundaries.
        assert_eq!(idx.overlap("chr1", 30, 35), (1, 6));
        assert_eq!(idx.overlap("chr1", 0, 5), (5, 6));
        // Query covering both spans plus the gap in between.
        assert_eq!(idx.overlap("chr1", 1, 30), (21, 30));
        // No overlap.
        assert_eq!(idx.overlap("chr1", 31, 40), (0, 10));
        assert_eq!(idx.overlap("chr1", 40, 50), (0, 11));
        // Chromosome absent from the index.
        assert_eq!(idx.overlap("chrX", 1, 10), (0, 10));
    }

    #[test]
    fn depth_at_least_basic() {
        let ivs = [(0u32, 10u32), (5, 15), (20, 25)];
        // depth 2 on [5,10), depth 1 on the rest of the union [0,15).
        assert_eq!(depth_at_least(&ivs, 1).to_string(), "0-14,20-24");
        assert_eq!(depth_at_least(&ivs, 2).to_string(), "5-9");
        assert_eq!(depth_at_least(&ivs, 3).to_string(), "-");
    }

    #[test]
    fn depth_at_least_adjacent_and_empty() {
        assert_eq!(depth_at_least(&[], 1).to_string(), "-");
        let adjacent = [(0u32, 10u32), (10, 20)];
        assert_eq!(depth_at_least(&adjacent, 1).to_string(), "0-19");
        let touching = [(0u32, 10u32), (10, 20), (10, 15)];
        assert_eq!(depth_at_least(&touching, 2).to_string(), "10-14");
    }

    #[test]
    fn depth_by_level_groups() {
        let ivs = [(0u32, 10u32), (5, 15)];
        let d = depth_by_level(&ivs, 1);
        assert_eq!(d["1"].to_string(), "0-4,10-14");
        assert_eq!(d["2"].to_string(), "5-9");
    }

    #[test]
    fn span_ops_fill_and_excise() {
        let mut s = BTreeMap::new();
        s.insert("chr1".to_string(), set("1-3,7-10,15-16"));
        // fill holes <= 3: 3-7 gap of 3 -> merged; 10-15 gap of 4 stays.
        assert_eq!(
            span_op(&s, SpanOp::Fill, 3)["chr1"].to_string(),
            "1-10,15-16"
        );
        // excise spans < 3: 15-16 (len 2) removed.
        assert_eq!(
            span_op(&s, SpanOp::Excise, 3)["chr1"].to_string(),
            "1-3,7-10"
        );
    }

    #[test]
    fn compare_intersect() {
        let mut first = BTreeMap::new();
        let mut s1 = BTreeMap::new();
        s1.insert("chr1".to_string(), set("1-10,20-30"));
        first.insert("a".to_string(), s1);
        let mut s2 = BTreeMap::new();
        s2.insert("chr1".to_string(), set("5-25"));
        let r = compare_sets(&first, &[s2], CompareOp::Intersect);
        assert_eq!(r["a"]["chr1"].to_string(), "5-10,20-25");
    }

    #[test]
    fn merge_keys_from_stems() {
        let dir = tempfile::tempdir().unwrap();
        let f1 = dir.path().join("sample.one.json");
        let f2 = dir.path().join("other.two.json");
        std::fs::write(&f1, r#"{"chr1":"1-5"}"#).unwrap();
        std::fs::write(&f2, r#"{"chr2":"6-9"}"#).unwrap();
        let m = merge_files(
            &[
                f1.to_string_lossy().into_owned(),
                f2.to_string_lossy().into_owned(),
            ],
            false,
        )
        .unwrap();
        assert_eq!(m.len(), 2);
        assert!(m.contains_key("sample"));
        assert!(m.contains_key("other"));
    }

    #[test]
    fn empty_json_is_single_empty_set() {
        let json = BTreeMap::new();
        let sets = json_to_sets(&json).unwrap();
        assert_eq!(sets.len(), 1);
        assert!(sets.contains_key("__single__"));
        assert!(sets["__single__"].is_empty());
    }

    #[test]
    fn json_to_set_rejects_invalid_values() {
        let mut json = BTreeMap::new();
        json.insert("chr1".to_string(), serde_json::Value::String("1-".into()));
        assert!(json_to_set(&json).is_err());
        let mut json = BTreeMap::new();
        json.insert(
            "chr1".to_string(),
            serde_json::Value::String("99999999999".into()),
        );
        assert!(json_to_set(&json).is_err());
        let mut json = BTreeMap::new();
        json.insert("chr1".to_string(), serde_json::json!({"x": "1-5"}));
        assert!(json_to_set(&json).is_err());
    }

    #[test]
    fn json_to_sets_rejects_mixed_shapes() {
        let mut json = BTreeMap::new();
        json.insert("a".to_string(), serde_json::json!({"chr1": "1-5"}));
        json.insert("b".to_string(), serde_json::Value::String("1-5".into()));
        // First value is an object (multi form): the string entry is invalid.
        assert!(json_to_sets(&json).is_err());

        let mut flat = BTreeMap::new();
        flat.insert("chr1".to_string(), serde_json::Value::String("1-5".into()));
        flat.insert("chr2".to_string(), serde_json::json!({"x": "1-5"}));
        // First value is a string (flat form): the object entry is invalid.
        assert!(json_to_sets(&flat).is_err());
    }

    #[test]
    fn combine_sets_unions_and_intersects() {
        let mut set_of = BTreeMap::new();
        let mut a = BTreeMap::new();
        a.insert("chr1".to_string(), set("1-10,20-30"));
        let mut b = BTreeMap::new();
        b.insert("chr1".to_string(), set("5-25"));
        b.insert("chr2".to_string(), set("1-50"));
        set_of.insert("a".to_string(), a);
        set_of.insert("b".to_string(), b);
        assert_eq!(
            combine_sets(&set_of, CompareOp::Union)["chr1"].to_string(),
            "1-30"
        );
        assert_eq!(
            combine_sets(&set_of, CompareOp::Intersect)["chr1"].to_string(),
            "5-10,20-25"
        );
        assert_eq!(
            combine_sets(&set_of, CompareOp::Intersect)["chr2"].to_string(),
            "-"
        );
    }

    #[test]
    fn convert_set_longest_and_all() {
        let mut runset = BTreeMap::new();
        runset.insert("chr1".to_string(), set("1-10,20-30"));
        assert_eq!(
            convert_set(&runset, false),
            vec!["chr1:1-10".to_string(), "chr1:20-30".to_string()]
        );
        assert_eq!(convert_set(&runset, true), vec!["chr1:20-30".to_string()]);
    }

    #[test]
    fn genome_set_covers_full_length() {
        let mut sizes = BTreeMap::new();
        sizes.insert("chr1".to_string(), 1000);
        sizes.insert("chr2".to_string(), 200);
        let g = genome_set(&sizes).unwrap();
        assert_eq!(g["chr1"].to_string(), "1-1000");
        assert_eq!(g["chr2"].to_string(), "1-200");
    }

    #[test]
    fn genome_set_rejects_non_positive_sizes() {
        let mut sizes = BTreeMap::new();
        sizes.insert("chr1".to_string(), 0);
        assert!(genome_set(&sizes).is_err());
        sizes.insert("chr2".to_string(), -5);
        assert!(genome_set(&sizes).is_err());
    }

    #[test]
    fn some_json_filters_keys() {
        let mut json = BTreeMap::new();
        json.insert("chr1".to_string(), serde_json::Value::String("1-5".into()));
        json.insert("chr2".to_string(), serde_json::Value::String("6-9".into()));
        let names: std::collections::BTreeSet<String> = ["chr1".to_string()].into_iter().collect();
        let out = some_json(&json, &names);
        assert_eq!(out.len(), 1);
        assert!(out.contains_key("chr1"));
    }

    #[test]
    fn split_json_requires_objects() {
        let mut json = BTreeMap::new();
        json.insert("a".to_string(), serde_json::json!({"chr1": "1-5"}));
        json.insert("b".to_string(), serde_json::Value::String("1-5".into()));
        assert!(split_json(&json).is_err());
        json.remove("b");
        let parts = split_json(&json).unwrap();
        assert_eq!(parts.len(), 1);
        assert_eq!(parts[0].0, "a");
    }

    #[test]
    fn stat_lines_and_statop() {
        let mut sizes = BTreeMap::new();
        sizes.insert("chr1".to_string(), 1000);
        let mut runset = BTreeMap::new();
        runset.insert("chr1".to_string(), set("1-500"));
        let s = stat_lines(&runset, &sizes, false, None).unwrap();
        assert_eq!(s, "chr1,1000,500,0.5000\nall,1000,500,0.5000");
        let s2 = runset.clone();
        let op = runset.clone();
        let so = statop_lines(&runset, &sizes, &s2, &op, false, None).unwrap();
        assert_eq!(
            so,
            "chr1,1000,500,500,500,0.5000,1.0000,2.0000\nall,1000,500,500,500,0.5000,1.0000,2.0000"
        );
        // A chromosome missing from sizes is a friendly error, not a panic.
        let mut bad = BTreeMap::new();
        bad.insert("chrX".to_string(), set("1-5"));
        assert!(stat_lines(&bad, &sizes, false, None).is_err());
    }
}
