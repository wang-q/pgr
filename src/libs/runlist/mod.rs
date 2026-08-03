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
use std::collections::BTreeMap;
use std::io::BufRead;

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
        if !range.is_valid() {
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
        if !range.is_valid() {
            continue;
        }
        iv_of
            .entry(range.chr().clone())
            .or_default()
            .push((*range.start() as u32, *range.end() as u32 + 1));
    }
    Ok(iv_of)
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
pub fn json_to_set(json: &BTreeMap<String, serde_json::Value>) -> BTreeMap<String, IntSpan> {
    json.iter()
        .filter_map(|(chr, v)| {
            let s = v.as_str()?;
            Some((chr.clone(), IntSpan::from(s)))
        })
        .collect()
}

/// Convert a runlist JSON map into `name -> chromosome -> IntSpan`; detects
/// the multi form (`{"name": {"chr": "..."}}`) automatically. An empty or
/// flat input is treated as a single set under the `__single__` key.
pub fn json_to_sets(
    json: &BTreeMap<String, serde_json::Value>,
) -> BTreeMap<String, BTreeMap<String, IntSpan>> {
    let is_multi = json.values().next().map(|v| v.is_object()).unwrap_or(false);
    if is_multi {
        json.iter()
            .filter_map(|(name, v)| {
                let inner = v.as_object()?;
                let set = json_to_set(&inner.iter().map(|(k, v)| (k.clone(), v.clone())).collect());
                Some((name.clone(), set))
            })
            .collect()
    } else {
        let mut m = BTreeMap::new();
        m.insert("__single__".to_string(), json_to_set(json));
        m
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
        let sets = json_to_sets(&json);
        assert_eq!(sets.len(), 1);
        assert!(sets.contains_key("__single__"));
        assert!(sets["__single__"].is_empty());
    }
}
