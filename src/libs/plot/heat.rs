//! GC-content × coverage heatmap rendering (KatGC heat plot equivalent).

use indexmap::IndexMap;
use std::io::Write;

/// Render-ready data for the LaTeX heatmap.
#[derive(Debug, Clone)]
pub struct GcHeatmap {
    /// `x y density` table rows (density = count / zmax, 0..1).
    pub table: String,
    /// X tick positions (count bins).
    pub xticks: Vec<f64>,
    /// X tick labels (coverage values).
    pub xtick_labels: Vec<String>,
    /// Y group labels (GC counts).
    pub ygroups: Vec<String>,
    /// Y tick positions.
    pub yticks: Vec<f64>,
    /// Figure width in cm.
    pub width: f64,
    /// Figure height in cm.
    pub height: f64,
    /// Longest y label width in ex.
    pub label_len: usize,
}

/// Build heatmap data from the raw matrix: 2x2 neighbor average, normalized
/// by the peak, with adaptive axis ticks.
pub fn heatmap(plot: &[Vec<u64>], xmax: usize, zmax: u64) -> GcHeatmap {
    let n_rows = plot.len().saturating_sub(1);
    let mut density = IndexMap::new();
    for i in 0..n_rows {
        let row0 = &plot[i];
        let row1 = &plot[i + 1];
        let values: Vec<f64> = (0..xmax)
            .map(|a| {
                let val = (row0[a] + row0[a + 1] + row1[a] + row1[a + 1]) / 4;
                (val.min(zmax) as f64) / (zmax.max(1) as f64)
            })
            .collect();
        density.insert(format!("GC {i}"), values);
    }
    density_to_heatmap(density, xmax)
}

/// Build heatmap data from `.kgc` rows `(gc, kf, count)` (already
/// 2x2-averaged and clamped by `gc`, which wrote the matrix).
pub fn heatmap_from_kgc(rows: &[(usize, usize, u64)], zmax: u64) -> GcHeatmap {
    let xmax = rows.iter().map(|&(_, a, _)| a).max().unwrap_or(0) + 1;
    let mut density = IndexMap::new();
    for (i, a, c) in rows {
        let row = density
            .entry(format!("GC {i}"))
            .or_insert_with(|| vec![0.0; xmax]);
        row[*a] = (*c as f64) / (zmax.max(1) as f64);
    }
    density_to_heatmap(density, xmax)
}

/// Shared axis construction from a per-GC-row density table.
fn density_to_heatmap(density: IndexMap<String, Vec<f64>>, xmax: usize) -> GcHeatmap {
    let table = super::histogram::create_table(&density);
    let n_groups = density.len().max(2);
    let ygroups: Vec<String> = density.keys().cloned().collect();
    let label_len = ygroups.iter().map(|s| s.len()).max().unwrap_or(3);

    // Adaptive ticks: ~8 ticks per axis, avoiding overcrowding.
    let xtick_step = (xmax.div_ceil(8)).max(1);
    let xticks: Vec<f64> = (0..xmax)
        .step_by(xtick_step)
        .map(|a| a as f64 - 0.5)
        .collect();
    let xtick_labels: Vec<String> = (0..xmax)
        .step_by(xtick_step)
        .map(|a| a.to_string())
        .collect();
    let ytick_step = (density.len().div_ceil(8)).max(1);
    let yticks: Vec<f64> = (0..=n_groups)
        .step_by(ytick_step)
        .map(|i| i as f64 - 0.5)
        .collect();

    // Cap the figure so wide matrices stay printable: cell 0.5 x 1.5 cm.
    let width = ((xmax as f64) * 0.5).min(25.0);
    let height = (n_groups as f64) * 1.5;
    GcHeatmap {
        table,
        xticks,
        xtick_labels,
        ygroups,
        yticks,
        width,
        height,
        label_len,
    }
}

/// Render a heatmap to a LaTeX writer (pgfplots, compile with tectonic).
pub fn render_heat<W: Write>(w: &mut W, hm: &GcHeatmap) -> anyhow::Result<()> {
    let mut context = tera::Context::new();
    context.insert("table", &hm.table);
    context.insert("xlabel", "k-mer coverage");
    context.insert("ylabel", "GC content");
    context.insert("width", &hm.width);
    context.insert("height", &hm.height);
    context.insert("xticks", &hm.xticks);
    context.insert("xtick_labels", &hm.xtick_labels);
    context.insert("ygroups", &hm.ygroups);
    context.insert("yticks", &hm.yticks);
    context.insert("label_len", &hm.label_len);
    super::histogram::render_hh_tex(&context, w)
}
