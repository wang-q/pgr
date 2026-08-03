//! Dot plot (collinear alignment plot) generation from PAF records.

use crate::libs::paf::PafRecord;
use std::collections::HashMap;

/// Options controlling dot plot rendering.
pub struct DotOpts {
    /// Minimum alignment block length to plot.
    pub min_len: u32,
    /// Minimum identity (`matches / block_length`) to plot.
    pub min_identity: f64,
    /// Identity at which the color scale saturates (deepest color).
    pub identity_max: f64,
    /// Maximum number of alignments to plot; `0` keeps all.
    pub max_align: usize,
    /// Plot area width in pixels; height is scaled from the axis lengths.
    pub width: u32,
    /// Make the plot frame square (independent x/y scaling).
    pub square: bool,
    /// Target-side region to zoom into (1-based inclusive); the query axis
    /// auto-focuses on the aligned regions.
    pub range: Option<PlotRange>,
}

/// A 1-based inclusive genomic region for axis clipping.
#[derive(Debug, Clone, PartialEq)]
pub struct PlotRange {
    /// Sequence name.
    pub chr: String,
    /// 1-based inclusive start.
    pub start: u32,
    /// 1-based inclusive end.
    pub end: u32,
}

impl std::str::FromStr for PlotRange {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> anyhow::Result<Self> {
        let r = crate::libs::ds::Range::from_str(s);
        if !r.is_valid() {
            anyhow::bail!("invalid range: {s}");
        }
        let start = *r.start() as u32;
        let end = *r.end() as u32;
        if start < 1 || end < start {
            anyhow::bail!("invalid range: {s}");
        }
        Ok(PlotRange {
            chr: r.chr().clone(),
            start,
            end,
        })
    }
}

/// Default `--min-len` value.
pub const DEFAULT_MIN_LEN: u32 = 100;
/// Default `--min-identity` value.
pub const DEFAULT_MIN_IDENTITY: f64 = 0.7;
/// Default `--max-identity` value.
pub const DEFAULT_MAX_IDENTITY: f64 = 1.0;
/// Default `--max-align` value.
pub const DEFAULT_MAX_ALIGN: usize = 100_000;
/// Default `--width` value.
pub const DEFAULT_WIDTH: u32 = 1200;

const TICK_TARGET_PX: f64 = 120.0;
/// Maximum gap (bp) between query alignments merged into one cluster.
const CLUSTER_GAP: i64 = 100_000;
/// A cluster is kept when its aligned bases are at least this fraction of the
/// largest cluster's aligned bases.
const CLUSTER_SCORE_RATIO: u64 = 100;

/// Per-cluster info used when rendering the query axis: display label (the
/// original sequence name) and true 0-based start of the cluster.
struct QueryClusterInfo {
    label: String,
    lo: u32,
}

/// Blues color ramp used for the identity color scale.
const BLUES: [(u8, u8, u8); 9] = [
    (247, 251, 255), // #F7FBFF
    (222, 235, 247), // #DEEBF7
    (198, 219, 239), // #C6DBEF
    (158, 202, 225), // #9ECAE1
    (107, 174, 214), // #6BAED6
    (66, 146, 198),  // #4292C6
    (33, 113, 181),  // #2171B5
    (8, 81, 156),    // #08519C
    (8, 48, 107),    // #08306B
];

/// Reds color ramp used for reverse-strand (non-collinear) alignments.
const REDS: [(u8, u8, u8); 9] = [
    (255, 245, 240), // #FFF5F0
    (254, 224, 210), // #FEE0D6
    (252, 187, 161), // #FCBBA1
    (252, 146, 114), // #FC9272
    (251, 106, 74),  // #FB6A4A
    (239, 59, 44),   // #EF3B2C
    (203, 24, 29),   // #CB181D
    (165, 15, 21),   // #A50F15
    (103, 0, 13),    // #67000D
];

/// Per-axis layout: contig order (first-appearance), lengths, offsets, and total length.
struct AxisInfo {
    names: Vec<String>,
    length: HashMap<String, u64>,
    offset: HashMap<String, u64>,
    total: u64,
}

fn build_axis<'a>(
    records: impl Iterator<Item = &'a PafRecord>,
    name_of: impl Fn(&PafRecord) -> &str,
    len_of: impl Fn(&PafRecord) -> u32,
) -> AxisInfo {
    let mut names: Vec<String> = Vec::new();
    let mut lens: HashMap<&str, u64> = HashMap::new();
    let mut order: HashMap<&str, usize> = HashMap::new();

    for rec in records {
        let name = name_of(rec);
        let len = len_of(rec) as u64;
        let e = lens.entry(name).or_insert(0);
        *e = (*e).max(len);
        if !order.contains_key(name) {
            order.insert(name, names.len());
            names.push(name.to_string());
        }
    }

    let mut length = HashMap::new();
    let mut offset = HashMap::new();
    let mut acc: u64 = 0;
    for name in &names {
        offset.insert(name.clone(), acc);
        let len = lens.get(name.as_str()).copied().unwrap_or(0);
        length.insert(name.clone(), len);
        acc += len;
    }

    AxisInfo {
        names,
        length,
        offset,
        total: acc,
    }
}

fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn identity(rec: &PafRecord) -> f64 {
    if rec.block_length == 0 {
        0.0
    } else {
        rec.matches as f64 / rec.block_length as f64
    }
}

/// "Nice" tick step (1/2/5 x 10^n) near `raw`.
fn nice_step(raw: f64) -> f64 {
    if raw <= 0.0 || !raw.is_finite() {
        return 1.0;
    }
    let mag = 10f64.powf(raw.log10().floor());
    let norm = raw / mag;
    let nice = if norm <= 1.0 {
        1.0
    } else if norm <= 2.0 {
        2.0
    } else if norm <= 5.0 {
        5.0
    } else {
        10.0
    };
    nice * mag
}

/// Compact bp label for tick marks (e.g. `1.5M`, `500k`, `1234`).
fn format_bp(v: u64) -> String {
    let k = (v as f64 / 1e3).round() as u64;
    if k >= 1000 {
        let m = format!("{:.2}", k as f64 / 1e3);
        let m = m.trim_end_matches('0').trim_end_matches('.');
        format!("{m}M")
    } else if v >= 1_000 {
        format!("{k}k")
    } else {
        v.to_string()
    }
}

/// Tick offsets (bp within a sequence) for an axis segment of `len` bp.
fn tick_positions(len: u64, scale: f64) -> Vec<u64> {
    let step = nice_step(TICK_TARGET_PX / scale.max(1e-9)).max(1.0) as u64;
    let mut ticks = Vec::new();
    let mut t = 0u64;
    while t < len {
        ticks.push(t);
        t = t.saturating_add(step);
    }
    ticks
}

/// Map `t` in [0, 1] to an SVG color on the given ramp.
fn ramp_color(t: f64, ramp: &[(u8, u8, u8); 9]) -> String {
    let t = t.clamp(0.0, 1.0);
    let scaled = t * (ramp.len() - 1) as f64;
    let i = (scaled.floor() as usize).min(ramp.len() - 2);
    let frac = scaled - i as f64;
    let (r1, g1, b1) = ramp[i];
    let (r2, g2, b2) = ramp[i + 1];
    let r = (r1 as f64 + (r2 as f64 - r1 as f64) * frac).round() as u8;
    let g = (g1 as f64 + (g2 as f64 - g1 as f64) * frac).round() as u8;
    let b = (b1 as f64 + (b2 as f64 - b1 as f64) * frac).round() as u8;
    format!("#{r:02X}{g:02X}{b:02X}")
}

/// Keep records overlapping `r` on the target axis, clipping their target
/// coordinates to `r` (0-based internally) and shrinking the axis to the
/// range length.
fn clip_target(records: Vec<PafRecord>, r: &PlotRange) -> Vec<PafRecord> {
    let pa_start = r.start - 1;
    let pa_end = r.end;
    let span = pa_end - pa_start;
    let mut out = Vec::new();
    for mut rec in records {
        if rec.target_name != r.chr {
            continue;
        }
        let s = rec.target_start.max(pa_start);
        let e = rec.target_end.min(pa_end);
        if s >= e {
            continue;
        }
        rec.target_start = s - pa_start;
        rec.target_end = e - pa_start;
        rec.target_length = span;
        out.push(rec);
    }
    out
}

/// Auto-focus the query axis on the significant aligned clusters: for each
/// query sequence, merge alignments greedily (gap <= `CLUSTER_GAP`) and keep
/// every cluster whose aligned bases are at least `1 / CLUSTER_SCORE_RATIO`
/// of the largest cluster (remote matches and matches on other chromosomes
/// stay visible; only tiny noise fragments are dropped). Each kept cluster
/// becomes its own axis segment, and its records are renamed to `chr#k` so
/// the shared axis layout treats them as separate segments. Rendering
/// subtracts the cluster start (negative offsets are clipped visually by the
/// SVG clipPath), so segment directions never distort.
fn focus_query(mut records: Vec<PafRecord>) -> (Vec<PafRecord>, HashMap<String, QueryClusterInfo>) {
    let mut by_chr: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, rec) in records.iter().enumerate() {
        by_chr.entry(rec.query_name.clone()).or_default().push(i);
    }

    let mut all_clusters: Vec<(String, Vec<usize>, u32, u32, u64)> = vec![];
    let mut max_score = 0u64;
    for (chr, mut idxs) in by_chr {
        idxs.sort_by_key(|&i| records[i].query_start);
        let mut cur = idxs[0];
        let mut members = vec![idxs[0]];
        let mut s = records[cur].query_start;
        let mut e = records[cur].query_end;
        let mut score = records[cur].block_length as u64;
        for &i in idxs.iter().skip(1) {
            if records[i].query_start as i64 - records[cur].query_end as i64 <= CLUSTER_GAP {
                e = e.max(records[i].query_end);
                score += records[i].block_length as u64;
                members.push(i);
                cur = i;
            } else {
                all_clusters.push((chr.clone(), members.clone(), s, e, score));
                max_score = max_score.max(score);
                members = vec![i];
                cur = i;
                s = records[i].query_start;
                e = records[i].query_end;
                score = records[i].block_length as u64;
            }
        }
        all_clusters.push((chr.clone(), members, s, e, score));
        max_score = max_score.max(score);
    }

    let mut info_of: HashMap<String, QueryClusterInfo> = HashMap::new();
    let mut seq_no: HashMap<String, usize> = HashMap::new();
    let mut keep = vec![true; records.len()];
    for (chr, members, s, e, score) in all_clusters {
        if score.saturating_mul(CLUSTER_SCORE_RATIO) < max_score {
            for i in &members {
                keep[*i] = false;
            }
            continue;
        }
        let k = seq_no.entry(chr.clone()).or_insert(0);
        let key = format!("{chr}#{k}");
        *k += 1;
        let span = e.saturating_sub(s).max(1);
        info_of.insert(
            key.clone(),
            QueryClusterInfo {
                label: chr.clone(),
                lo: s,
            },
        );
        for i in members {
            records[i].query_name = key.clone();
            records[i].query_length = span;
        }
    }
    records = records
        .into_iter()
        .zip(keep)
        .filter_map(|(r, k)| k.then_some(r))
        .collect();
    (records, info_of)
}

/// Render a dot plot as an SVG string from PAF records.
pub fn render_dot_svg(records: &[PafRecord], opts: &DotOpts) -> anyhow::Result<String> {
    if opts.identity_max <= opts.min_identity {
        anyhow::bail!(
            "identity_max ({}) must be greater than min_identity ({})",
            opts.identity_max,
            opts.min_identity
        );
    }

    let mut filtered: Vec<PafRecord> = records
        .iter()
        .filter(|r| r.block_length >= opts.min_len && identity(r) >= opts.min_identity)
        .cloned()
        .collect();

    if filtered.is_empty() {
        anyhow::bail!(
            "no alignments pass filters (min_len={}, min_identity={})",
            opts.min_len,
            opts.min_identity
        );
    }

    if opts.max_align > 0 && filtered.len() > opts.max_align {
        let mut idx: Vec<usize> = (0..filtered.len()).collect();
        idx.sort_by(|&a, &b| filtered[b].block_length.cmp(&filtered[a].block_length));
        idx.truncate(opts.max_align);
        let mut keep = vec![false; filtered.len()];
        for i in idx {
            keep[i] = true;
        }
        filtered = filtered
            .into_iter()
            .zip(keep)
            .filter_map(|(r, k)| k.then_some(r))
            .collect();
    }

    // True 0-based coordinate of the start of each axis segment (zoom offsets).
    let x_offset: u64 = opts.range.as_ref().map(|r| r.start as u64 - 1).unwrap_or(0);
    let mut q_offset_of: HashMap<String, QueryClusterInfo> = HashMap::new();
    if let Some(r) = &opts.range {
        filtered = clip_target(filtered, r);
        if filtered.is_empty() {
            anyhow::bail!("no alignments intersect the requested range");
        }
        let (focused, lo_of) = focus_query(filtered);
        filtered = focused;
        q_offset_of = lo_of;
    }

    let t_axis = build_axis(filtered.iter(), |r| &r.target_name, |r| r.target_length);
    let q_axis = build_axis(filtered.iter(), |r| &r.query_name, |r| r.query_length);

    if t_axis.total == 0 || q_axis.total == 0 {
        anyhow::bail!("alignments have zero-length axis sequences");
    }

    let plot_w = opts.width.max(1) as f64;
    let x_scale = plot_w / t_axis.total as f64;
    let y_scale = if opts.square {
        plot_w / q_axis.total as f64
    } else {
        x_scale
    };
    let plot_h = if opts.square {
        plot_w
    } else {
        q_axis.total as f64 * x_scale
    };
    let line_width = (opts.width.max(1) as f64 / 200.0).clamp(1.0, 20.0);
    let label_font = (plot_w / 60.0).clamp(12.0, 32.0);
    let tick_font = label_font * 0.8;
    // Short segments are extended to at least 1.5x the stroke width so they
    // stay line-like instead of collapsing into dots narrower than the line.
    let min_line_px = (line_width * 1.5).max(4.0);
    let mut margin_top = (plot_w * 0.02).clamp(12.0, 40.0);

    // If the topmost query segment is too short for its rotated name, grow
    // the top margin so the name stays inside the viewBox.
    if let Some(last) = q_axis.names.last() {
        let seg_px = q_axis.length[last] as f64 * y_scale;
        let label = q_offset_of
            .get(last)
            .map(|i| i.label.chars().count())
            .unwrap_or_else(|| last.chars().count());
        let name_px = label as f64 * label_font;
        if seg_px < name_px {
            margin_top = margin_top.max(name_px / 2.0 - seg_px / 2.0 + 12.0);
        }
    }
    // Keep room for the legend above the plot frame.
    margin_top = margin_top.max(90.0);
    let margin_right = (plot_w * 0.02).clamp(12.0, 40.0);
    let plot_top = margin_top;
    let plot_bottom = plot_top + plot_h;

    // Y-axis names are rotated -90 degrees; overlapping ones move to columns
    // to the left, with alternating black/dark-gray colors.
    let mut y_place: Vec<(String, f64, f64, usize)> = Vec::new();
    let mut y_cols: Vec<Vec<(f64, f64)>> = Vec::new();
    for name in &q_axis.names {
        let off = q_axis.offset[name];
        let len = q_axis.length[name];
        if len as f64 * y_scale < label_font {
            continue;
        }
        let cy = plot_bottom - (off as f64 + len as f64 / 2.0) * y_scale;
        let label = q_offset_of
            .get(name)
            .map(|i| i.label.as_str())
            .unwrap_or(name);
        let name_px = label.chars().count() as f64 * label_font;
        let mut col = 0;
        loop {
            if col == y_cols.len() {
                y_cols.push(Vec::new());
            }
            let overlap = y_cols[col]
                .iter()
                .any(|&(a, b)| cy - name_px / 2.0 < b && cy + name_px / 2.0 > a);
            if !overlap {
                y_cols[col].push((cy - name_px / 2.0, cy + name_px / 2.0));
                break;
            }
            col += 1;
        }
        y_place.push((label.to_string(), cy, name_px, col));
    }

    let mut margin_left = tick_font * 3.0 + label_font + 28.0;
    if y_cols.len() > 1 {
        margin_left = margin_left.max(
            8.0 + tick_font * 3.0
                + label_font
                + 10.0
                + (y_cols.len() - 1) as f64 * (label_font + 8.0)
                + 6.0,
        );
    }
    let plot_left = margin_left;
    let name_x = plot_left - 8.0 - tick_font * 3.0 - label_font / 2.0 - 10.0;

    // X-axis names are horizontal, alternating black/dark-gray per segment.
    let mut x_place: Vec<(String, f64)> = Vec::new();
    for name in &t_axis.names {
        let off = t_axis.offset[name];
        let len = t_axis.length[name];
        if len as f64 * x_scale < label_font {
            continue;
        }
        let cx = plot_left + (off as f64 + len as f64 / 2.0) * x_scale;
        x_place.push((name.clone(), cx));
    }

    let flat_y = plot_bottom + tick_font + 8.0 + tick_font + 6.0 + label_font;
    let mut margin_bottom = tick_font * 3.0 + label_font * 2.0 + 50.0;
    margin_bottom = margin_bottom.max(flat_y - plot_bottom + label_font + 10.0);
    let svg_w = margin_left + plot_w + margin_right;
    let svg_h = plot_top + plot_h + margin_bottom;

    let mut segments = String::new();
    for rec in &filtered {
        let mut x0 = plot_left
            + (t_axis.offset[&rec.target_name] + rec.target_start as u64) as f64 * x_scale;
        let mut x1 =
            plot_left + (t_axis.offset[&rec.target_name] + rec.target_end as u64) as f64 * x_scale;
        // PAF stores reverse-strand query coordinates in forward orientation,
        // so swap them to draw the segment with a negative slope.
        let (qa, qb) = if rec.strand == '-' {
            (rec.query_end, rec.query_start)
        } else {
            (rec.query_start, rec.query_end)
        };
        // Flip the y axis so query position 0 sits at the bottom and forward
        // alignments run from bottom-left to top-right.
        let mut y0 = plot_top + plot_h
            - (q_axis.offset[&rec.query_name] as i64 + qa as i64
                - q_offset_of
                    .get(&rec.query_name)
                    .map(|i| i.lo as i64)
                    .unwrap_or(0)) as f64
                * y_scale;
        let mut y1 = plot_top + plot_h
            - (q_axis.offset[&rec.query_name] as i64 + qb as i64
                - q_offset_of
                    .get(&rec.query_name)
                    .map(|i| i.lo as i64)
                    .unwrap_or(0)) as f64
                * y_scale;
        // Sub-pixel segments are invisible; extend them to a minimum visible length.
        let dx = x1 - x0;
        let dy = y1 - y0;
        let len = (dx * dx + dy * dy).sqrt();
        if len > 0.0 && len < min_line_px {
            let ext = (min_line_px - len) / 2.0;
            let (ux, uy) = (dx / len, dy / len);
            x0 -= ux * ext;
            x1 += ux * ext;
            y0 -= uy * ext;
            y1 += uy * ext;
        }
        let t = (identity(rec) - opts.min_identity) / (opts.identity_max - opts.min_identity);
        let color = if rec.strand == '+' {
            ramp_color(t, &BLUES)
        } else {
            ramp_color(t, &REDS)
        };
        let line = format!(
            r#"<line x1="{x0:.1}" y1="{y0:.1}" x2="{x1:.1}" y2="{y1:.1}" stroke="{color}"/>"#
        );
        segments.push_str(&line);
    }

    // Grid lines and tick marks for both axes, plus contig separators.
    let mut grid = String::new();
    let mut separators = String::new();
    let mut x_ticks = String::new();
    let mut y_ticks = String::new();

    for name in &t_axis.names {
        let off = t_axis.offset[name];
        let len = t_axis.length[name];
        let x0 = plot_left + off as f64 * x_scale;
        separators.push_str(&format!(
            r#"<line x1="{x0:.1}" y1="{plot_top:.1}" x2="{x0:.1}" y2="{plot_bottom:.1}"/>"#
        ));
        // Segments too short for tick labels keep only the separator line.
        if len as f64 * x_scale >= tick_font {
            for t in tick_positions(len, x_scale) {
                let x = plot_left + (off as f64 + t as f64) * x_scale;
                grid.push_str(&format!(
                    r#"<line x1="{x:.1}" y1="{plot_top:.1}" x2="{x:.1}" y2="{plot_bottom:.1}"/>"#
                ));
                // Each axis segment starts its own scale at 0; in zoom mode
                // the single segment shows true genomic coordinates.
                let label = if opts.range.is_some() {
                    format_bp(x_offset + off + t)
                } else {
                    format_bp(t)
                };
                x_ticks.push_str(&format!(
                    r#"<text x="{x:.1}" y="{}" font-size="{tick_font:.0}" text-anchor="middle">{}</text>"#,
                    plot_bottom + tick_font + 8.0,
                    label
                ));
            }
        }
    }

    for name in &q_axis.names {
        let off = q_axis.offset[name];
        let len = q_axis.length[name];
        let q_lo = q_offset_of.get(name).map(|i| i.lo as u64).unwrap_or(0);
        let y0 = plot_bottom - off as f64 * y_scale;
        separators.push_str(&format!(
            r#"<line x1="{plot_left:.1}" y1="{y0:.1}" x2="{:.1}" y2="{y0:.1}"/>"#,
            plot_left + plot_w
        ));
        // Segments too short for tick labels keep only the separator line.
        if len as f64 * y_scale >= tick_font {
            for t in tick_positions(len, y_scale) {
                let y = plot_bottom - (off as f64 + t as f64) * y_scale;
                grid.push_str(&format!(
                    r#"<line x1="{plot_left:.1}" y1="{y:.1}" x2="{:.1}" y2="{y:.1}"/>"#,
                    plot_left + plot_w
                ));
                let label = format_bp(q_lo + t);
                y_ticks.push_str(&format!(
                    r#"<text x="{:.1}" y="{y:.1}" font-size="{tick_font:.0}" text-anchor="end">{}</text>"#,
                    plot_left - 8.0,
                    label
                ));
            }
        }
    }

    let mut x_labels = String::new();
    for (i, (name, cx)) in x_place.iter().enumerate() {
        let color = if i % 2 == 0 { "#000000" } else { "#777777" };
        x_labels.push_str(&format!(
            r##"<text x="{cx:.1}" y="{flat_y:.1}" font-size="{label_font:.0}" font-weight="bold" fill="{color}" text-anchor="middle">{}</text>"##,
            escape_xml(name)
        ));
    }

    let mut y_labels = String::new();
    for (i, (label, cy, _, col)) in y_place.iter().enumerate() {
        let color = if i % 2 == 0 { "#000000" } else { "#777777" };
        let x = name_x - *col as f64 * (label_font + 8.0);
        y_labels.push_str(&format!(
            r##"<text x="{x:.1}" y="{cy:.1}" font-size="{label_font:.0}" font-weight="bold" fill="{color}" text-anchor="middle" transform="rotate(-90 {x:.1} {cy:.1})">{}</text>"##,
            escape_xml(label)
        ));
    }

    let mut out = String::new();
    out.push_str(&format!(
        r#"<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" font-family="sans-serif" width="{svg_w:.0}" height="{svg_h:.0}" viewBox="0 0 {svg_w:.0} {svg_h:.0}">
"#
    ));
    out.push_str(r#"<rect width="100%" height="100%" fill="white"/>"#);
    out.push_str(&format!(
        r##"<defs><clipPath id="plot-clip"><rect x="{plot_left:.1}" y="{plot_top:.1}" width="{plot_w:.1}" height="{plot_h:.1}"/></clipPath></defs>"##
    ));
    out.push_str(&format!(
        r##"<g id="grid" stroke="#EDEDED" stroke-width="{:.2}">{grid}</g>"##,
        line_width * 0.3
    ));
    out.push_str(&format!(
        r##"<rect x="{:.1}" y="{:.1}" width="{plot_w:.1}" height="{plot_h:.1}" fill="none" stroke="#808080" stroke-width="{line_width:.1}"/>"##,
        plot_left,
        plot_top
    ));
    out.push_str(&format!(
        r##"<g id="separators" stroke="#808080" stroke-width="{line_width:.1}">{separators}</g>"##
    ));
    out.push_str(&format!(
        r##"<g id="segments" stroke-width="{line_width:.1}" opacity="0.8" clip-path="url(#plot-clip)">{segments}</g>"##
    ));
    out.push_str(&x_ticks);
    out.push_str(&y_ticks);
    out.push_str(&x_labels);
    out.push_str(&y_labels);
    out.push_str(&legend(plot_left, plot_top, plot_w, label_font, opts));
    out.push_str("</svg>\n");

    Ok(out)
}

/// Build the identity color-scale legend above the plot's top-right corner.
fn legend(plot_left: f64, plot_top: f64, plot_w: f64, label_font: f64, opts: &DotOpts) -> String {
    let legend_w = 260.0;
    let legend_h = 14.0;
    // Shift the legend left so the 100% label stays inside the plot frame.
    let lx = plot_left + plot_w - legend_w - 80.0;
    let ly = plot_top - 46.0;
    let label_y = ly + legend_h / 2.0 + label_font * 0.5;
    let lo = format!("{:.0}%", opts.min_identity * 100.0);
    let hi = format!("{:.0}%", opts.identity_max * 100.0);

    let mut s = String::new();
    s.push_str(&format!(
        r#"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="62" rx="4" fill="white" opacity="0.85"/>"#,
        lx - 44.0,
        plot_top - 86.0,
        legend_w + 96.0
    ));
    const STEPS: usize = 12;
    for i in 0..STEPS {
        let t = i as f64 / (STEPS - 1) as f64;
        let color = ramp_color(t, &BLUES);
        let x = lx + i as f64 * legend_w / STEPS as f64;
        s.push_str(&format!(
            r#"<rect x="{x:.1}" y="{ly:.1}" width="{:.1}" height="{legend_h:.1}" fill="{color}"/>"#,
            legend_w / STEPS as f64
        ));
    }
    s.push_str(&format!(
        r#"<text x="{:.1}" y="{}" font-size="{:.0}" text-anchor="middle">Forward blue / Reverse red</text>"#,
        lx + legend_w / 2.0,
        plot_top - 62.0,
        label_font
    ));
    s.push_str(&format!(
        r#"<text x="{:.1}" y="{label_y:.1}" font-size="{label_font:.0}" text-anchor="end">{lo}</text>"#,
        lx - 8.0
    ));
    s.push_str(&format!(
        r#"<text x="{:.1}" y="{label_y:.1}" font-size="{label_font:.0}">{hi}</text>"#,
        lx + legend_w + 8.0
    ));
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[allow(clippy::too_many_arguments)]
    fn rec(
        qname: &str,
        qstart: u32,
        qend: u32,
        strand: char,
        tname: &str,
        tstart: u32,
        tend: u32,
        matches: u32,
        block: u32,
    ) -> PafRecord {
        PafRecord {
            query_name: qname.to_string(),
            query_length: 1000,
            query_start: qstart,
            query_end: qend,
            strand,
            target_name: tname.to_string(),
            target_length: 1000,
            target_start: tstart,
            target_end: tend,
            matches,
            block_length: block,
            mapq: 60,
            tags: vec![],
        }
    }

    #[test]
    fn renders_segments_with_identity_colors() {
        let records = vec![
            rec("q1", 10, 500, '+', "t1", 20, 480, 450, 490),
            rec("q1", 600, 900, '-', "t2", 100, 350, 280, 300),
        ];
        let svg = render_dot_svg(
            &records,
            &DotOpts {
                min_len: 100,
                min_identity: 0.7,
                identity_max: 1.0,
                max_align: 0,
                width: 2000,
                square: false,
                range: None,
            },
        )
        .unwrap();
        assert!(svg.contains(r#"id="segments""#));
        let seg = svg.split(r#"id="segments""#).nth(1).unwrap();
        assert_eq!(seg.matches("<line ").count(), 2);
        // 0.918 (450/490) and 0.933 (280/300) map to different blues
        let colors: Vec<&str> = svg
            .match_indices("stroke=\"#")
            .map(|(i, _)| &svg[i + 8..i + 14])
            .collect();
        assert_ne!(colors[0], colors[1]);
        // legend labels use the width-scaled font (2000 / 60 = 33.3, clamped to 32)
        assert!(svg.contains("font-size=\"32\""));
    }

    #[test]
    fn filters_by_len_identity_and_max() {
        let records = vec![
            rec("q1", 0, 200, '+', "t1", 0, 200, 200, 200), // identity 1.0, len 200
            rec("q2", 0, 50, '+', "t2", 0, 50, 50, 50),     // len < min_len
            rec("q3", 0, 300, '+', "t3", 0, 300, 150, 300), // identity 0.5 < 0.7
            rec("q4", 0, 400, '-', "t4", 0, 400, 400, 400), // longest
        ];
        let svg = render_dot_svg(
            &records,
            &DotOpts {
                min_len: 100,
                min_identity: 0.7,
                identity_max: 1.0,
                max_align: 1,
                width: 500,
                square: false,
                range: None,
            },
        )
        .unwrap();
        // only the longest record (q4/t4, '-' strand) survives
        let seg = svg.split(r#"id="segments""#).nth(1).unwrap();
        assert_eq!(seg.matches("<line ").count(), 1);
        assert!(svg.contains("stroke=\"#"));
    }

    #[test]
    fn rejects_identity_max_le_min() {
        let records = vec![rec("q1", 0, 200, '+', "t1", 0, 200, 200, 200)];
        let err = render_dot_svg(
            &records,
            &DotOpts {
                min_len: 100,
                min_identity: 0.9,
                identity_max: 0.9,
                max_align: 0,
                width: 500,
                square: false,
                range: None,
            },
        )
        .unwrap_err();
        assert!(err.to_string().contains("must be greater"));
    }

    #[test]
    fn zooms_with_range_and_focuses_query() {
        let records = vec![
            rec("q1", 0, 1000, '+', "t1", 0, 1000, 1000, 1000),
            rec("q2", 200, 600, '-', "t1", 500, 900, 400, 400),
            rec("q3", 0, 800, '+', "t2", 0, 800, 800, 800),
        ];
        let range = "t1:100-600".parse::<PlotRange>().unwrap();
        let svg = render_dot_svg(
            &records,
            &DotOpts {
                min_len: 100,
                min_identity: 0.7,
                identity_max: 1.0,
                max_align: 0,
                width: 500,
                square: false,
                range: Some(range),
            },
        )
        .unwrap();
        // t2 record dropped; two records clipped into the range remain
        let seg = svg.split(r#"id="segments""#).nth(1).unwrap();
        assert_eq!(seg.matches("<line ").count(), 2);
        assert!(svg.contains(">t1</text>"));
        assert!(!svg.contains(">t2</text>"));
        assert!(svg.contains(">q1</text>"));
        assert!(svg.contains(">q2</text>"));
        assert!(!svg.contains(">q3</text>"));
    }

    #[test]
    fn square_makes_plot_frame_square() {
        let records = vec![
            rec("q1", 0, 1000, '+', "t1", 0, 1000, 1000, 1000),
            rec("q1", 1000, 2000, '+', "t1", 1000, 2000, 1000, 1000),
        ];
        let svg = render_dot_svg(
            &records,
            &DotOpts {
                min_len: 100,
                min_identity: 0.7,
                identity_max: 1.0,
                max_align: 0,
                width: 800,
                square: true,
                range: None,
            },
        )
        .unwrap();
        // the plot-frame rect must be square: width == height
        let re = regex::Regex::new(
            r#"<rect x="[0-9.]+" y="[0-9.]+" width="([0-9.]+)" height="([0-9.]+)" fill="none""#,
        )
        .unwrap();
        let caps = re.captures(&svg).expect("missing plot frame rect");
        let w: f64 = caps[1].parse().unwrap();
        let h: f64 = caps[2].parse().unwrap();
        assert!((w - h).abs() < 0.001, "frame {w} x {h}");
    }

    #[test]
    fn strand_determines_segment_slope() {
        // Forward strand: bottom-left to top-right; reverse: top-left to bottom-right.
        let records = vec![
            rec("q1", 100, 900, '+', "t1", 100, 900, 800, 800),
            rec("q2", 100, 900, '-', "t2", 100, 900, 800, 800),
        ];
        let svg = render_dot_svg(
            &records,
            &DotOpts {
                min_len: 100,
                min_identity: 0.7,
                identity_max: 1.0,
                max_align: 0,
                width: 1200,
                square: false,
                range: None,
            },
        )
        .unwrap();
        let re = regex::Regex::new(
            r##"<line x1="([0-9.]+)" y1="([0-9.]+)" x2="([0-9.]+)" y2="([0-9.]+)" stroke="#([0-9A-F]{6})""##,
        )
        .unwrap();
        let mut fwd_slope: Option<bool> = None;
        let mut rev_slope: Option<bool> = None;
        for c in re.captures_iter(&svg) {
            let y1: f64 = c[2].parse().unwrap();
            let y2: f64 = c[4].parse().unwrap();
            let color = &c[5];
            let is_blue = u8::from_str_radix(&color[4..6], 16).unwrap()
                > u8::from_str_radix(&color[0..2], 16).unwrap();
            let negative = y2 > y1;
            if is_blue {
                fwd_slope = Some(negative);
            } else {
                rev_slope = Some(negative);
            }
        }
        assert_eq!(fwd_slope, Some(false), "forward must go up-right");
        assert_eq!(rev_slope, Some(true), "reverse must go down-right");
    }

    #[test]
    fn errors_when_nothing_passes() {
        let records = vec![rec("q1", 0, 50, '+', "t1", 0, 50, 50, 50)];
        let err = render_dot_svg(
            &records,
            &DotOpts {
                min_len: 100,
                min_identity: 0.7,
                identity_max: 1.0,
                max_align: 0,
                width: 500,
                square: false,
                range: None,
            },
        )
        .unwrap_err();
        assert!(err.to_string().contains("no alignments pass"));
    }
}
