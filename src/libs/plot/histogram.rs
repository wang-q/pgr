use indexmap::IndexMap;

/// Axis-layout parameters for a histo-heatmap (hh) plot.
pub struct HhAxis {
    pub width: f64,
    pub height: f64,
    pub ygroups: Vec<String>,
    pub yticks: Vec<f64>,
    pub label_len: usize,
    pub xticks: Vec<f64>,
    pub xtick_labels: Vec<String>,
}

/// Compute the hh axis layout (width/height/ticks/labels) from density data, bin edges, and cell unit.
pub fn compute_hh_axis(
    density_data: &IndexMap<String, Vec<f64>>,
    bins: usize,
    bin_edges: &[f64],
    unit: (f64, f64),
) -> HhAxis {
    let width = (bins as f64) * unit.0;
    let height = (density_data.len() as f64).max(2.0) * unit.1;

    let ygroups: Vec<String> = density_data.keys().cloned().collect();
    let yticks = (0..=density_data.len().max(2))
        .map(|i| i as f64 - 0.5)
        .collect::<Vec<_>>();
    let label_len = ygroups.iter().map(|s| s.len()).max().unwrap_or(0).max(3);

    let xticks = (0..=bins)
        .filter_map(|i| {
            if i % 5 == 0 {
                Some(i as f64 - 0.5)
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    let xtick_labels = bin_edges
        .iter()
        .enumerate()
        .filter_map(|(i, &edge)| {
            if i % 5 == 0 {
                Some(format!("{}", edge))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();

    HhAxis {
        width,
        height,
        ygroups,
        yticks,
        label_len,
        xticks,
        xtick_labels,
    }
}

/// Render the hh LaTeX document from `context` and write it to `writer`.
pub fn render_hh_tex<W: std::io::Write>(
    context: &tera::Context,
    writer: &mut W,
) -> anyhow::Result<()> {
    static FILE_TEMPLATE: &str = include_str!("../../assets/heatmap.tex");
    let mut template = FILE_TEMPLATE.to_string();

    let out_string = r###"%
width={{ width }}cm,
height={{ height }}cm,
xlabel={ {{ xlabel }} },
ylabel={ {{ ylabel }} },
extra x ticks={ {{ xticks | join(sep=", ") }} },
extra x tick labels={ {{ xtick_labels | join(sep=", ") }} },
yticklabels={ {{ ygroups | join(sep=", ") }} },
extra y ticks={ {{ yticks | join(sep=", ") }} },
y tick label style={
    text width={{ label_len }}ex,
},
    "###;
    crate::libs::plot::common::replace_section(
        &mut template,
        "%AXIS_BEGIN",
        "%AXIS_END",
        out_string,
    )?;

    let table = crate::libs::plot::common::context_get_str(context, "table")?;
    crate::libs::plot::common::replace_section(&mut template, "%TABLE_BEGIN", "%TABLE_END", table)?;

    crate::libs::plot::common::render_and_write(&template, context, writer)?;
    Ok(())
}

/// Load grouped numeric data from a TSV file.
///
/// Reads column `col` (1-based) as numeric values, optionally grouping by
/// the column `group` (1-based). Returns `(data, xlabel, ylabel)` where
/// `data` maps group name → vector of values, `xlabel` is the header of the
/// value column, and `ylabel` is the header of the group column (empty if
/// no grouping).
pub fn load_data(
    infile: &str,
    col: usize,
    group: Option<&usize>,
) -> anyhow::Result<(IndexMap<String, Vec<f64>>, String, String)> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(crate::libs::io::reader(infile)?);

    let headers = rdr.headers()?.clone();
    let mut data: IndexMap<String, Vec<f64>> = IndexMap::new();

    let xlabel = headers[col - 1].to_string();
    let ylabel = match group {
        Some(g) => headers[*g - 1].to_string(),
        None => String::new(),
    };

    for result in rdr.records() {
        let record = result?;

        if let Ok(val) = record[col - 1].parse::<f64>() {
            let group_name = match group {
                Some(g) => record[*g - 1].to_string(),
                None => "default".to_string(),
            };

            data.entry(group_name).or_default().push(val);
        }
    }

    Ok((data, xlabel, ylabel))
}

/// Compute per-group histograms and bin edges from grouped values.
#[allow(clippy::type_complexity)]
pub fn calc_hist(
    data: &IndexMap<String, Vec<f64>>,
    bins: usize,
    xmm: Option<(f64, f64)>,
) -> anyhow::Result<(IndexMap<String, Vec<usize>>, Vec<f64>)> {
    // Calculate global range
    let (min_val, max_val) = match xmm {
        Some((min, max)) => (min, max),
        None => {
            let (min, max) = data
                .values()
                .flatten()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), &val| {
                    (min.min(val), max.max(val))
                });
            // Normalize range to neat values
            let magnitude_min = if min.abs() < f64::EPSILON {
                1.0
            } else {
                10f64.powf(min.abs().log10().floor())
            };
            let magnitude_max = if max.abs() < f64::EPSILON {
                1.0
            } else {
                10f64.powf(max.abs().log10().floor())
            };

            let norm_min = (min / magnitude_min).floor() * magnitude_min;
            let norm_max = (max / magnitude_max).ceil() * magnitude_max;

            (norm_min, norm_max)
        }
    };

    let mut hist_data = IndexMap::new();
    let bin_width = (max_val - min_val) / (bins as f64);

    // Calculate histogram for each group
    for (group_name, values) in data.iter() {
        let mut hist = vec![0usize; bins];
        for &val in values {
            if val >= min_val && val <= max_val {
                let bin = ((val - min_val) / bin_width).floor() as usize;
                let bin = bin.min(bins - 1); // Handle edge case
                hist[bin] += 1;
            }
        }
        hist_data.insert(group_name.clone(), hist);
    }

    // Calculate bin edges
    let mut bin_edges = Vec::with_capacity(bins + 1);
    for i in 0..=bins {
        bin_edges.push(min_val + (i as f64) * bin_width);
    }

    Ok((hist_data, bin_edges))
}

/// Convert per-group histogram counts into densities (fractions of group totals).
pub fn calc_density(hist_data: &IndexMap<String, Vec<usize>>) -> IndexMap<String, Vec<f64>> {
    let mut density_data = IndexMap::new();

    for (group_name, hist) in hist_data.iter() {
        let total_samples = hist.iter().sum::<usize>() as f64;
        let density: Vec<f64> = hist
            .iter()
            .map(|&count| (count as f64) / total_samples)
            .collect();
        density_data.insert(group_name.clone(), density);
    }

    density_data
}

/// Render density data as a TSV table of `x y density` rows for LaTeX heatmap.
pub fn create_table(density_data: &IndexMap<String, Vec<f64>>) -> String {
    let mut table = String::new();
    let bins = density_data.values().next().map_or(0, |v| v.len());

    // Iterate through each group
    for (y, (_, densities)) in density_data.iter().enumerate() {
        // Iterate through each bin
        for (x, &d) in densities.iter().enumerate().take(bins) {
            table.push_str(&format!(
                "    {:3} {:3} {:.4}\n",
                x, // x coordinate (3 digits)
                y, // y coordinate (3 digits)
                d  // density value (4 decimal places)
            ));
        }
        table.push('\n');
    }

    // Add a dummy group with zeros if there's only one group
    if density_data.len() == 1 {
        for x in 0..bins {
            table.push_str(&format!(
                "    {:3} {:3} {:.4}\n",
                x,   // x coordinate
                1,   // y coordinate (second group)
                0.0  // density value (zero)
            ));
        }
        table.push('\n');
    }

    table
}
