//! GenomeScope-style k-mer spectra rendering (genescopefk.R plots).
//!
//! Renders the four views (linear / transformed / log / transformed-log)
//! as a single standalone LaTeX document: observed histogram, full model,
//! error region, and k-mer peak markers, with the model summary in the
//! figure title.

use crate::libs::kmer::genomescope::{predict, ModelParams};
use std::io::Write;

/// Model summary shown in the figure title.
#[derive(Debug, Clone, Copy)]
pub struct SpectraSummary {
    /// K-mer length.
    pub k: usize,
    /// Ploidy.
    pub p: usize,
    /// Fitted parameters.
    pub params: ModelParams,
    /// Kmercov standard error (for the error-region cutoff).
    pub se_kmercov: f64,
    /// Genome haploid length (bp).
    pub len: f64,
    /// Unique-sequence fraction of the haploid length (R title
    /// `unique_len/total_len`; 1 when repeats are absent).
    pub unique: f64,
    /// Read error rate.
    pub error_rate: f64,
}

/// Sequencing error rate, mirroring genescopefk.R: the fraction of k-mer
/// instances unexplained by the model in the low-coverage error region,
/// converted to a per-base rate via `1 - (1 - ratio)^(1/k)`.
pub fn compute_error_rate(
    hist: &[(f64, f64)],
    k: usize,
    params: &ModelParams,
    se_kmercov: f64,
) -> f64 {
    let kcovfloor = ((params.kmercov - 2.0 * se_kmercov).floor().max(1.0)) as usize;
    // R `error_xcutoff_ind = tail(which(x <= error_xcutoff), 1)`.
    let err_cut = hist
        .iter()
        .rposition(|&(x, _)| x <= kcovfloor as f64)
        .map_or(1, |i| i + 1);
    let mut error_kmers = vec![0.0f64; err_cut];
    let mut first_zero = err_cut;
    for i in 0..err_cut {
        let (x, count) = hist[i];
        let pred = params.length * predict(1, 1, params, k, x);
        let e = count - pred;
        if e < 1.0 {
            first_zero = i;
            break;
        }
        error_kmers[i] = e;
    }
    for v in error_kmers[..first_zero].iter_mut() {
        *v = v.max(1e-10);
    }
    let total_error: f64 = error_kmers[..first_zero]
        .iter()
        .zip(&hist[..first_zero])
        .map(|(&e, &(x, _))| e * x)
        .sum();
    let total: f64 = hist.iter().map(|&(x, c)| x * c).sum();
    1.0 - (1.0 - total_error / total.max(1.0)).powf(1.0 / k as f64)
}

/// Render the four-view spectra document to `w`.
///
/// `hist` holds `(coverage, count)` pairs (zero-count rows absent, matching
/// the reference histogram file). The document compiles with tectonic.
pub fn render_spectra<W: Write>(
    w: &mut W,
    hist: &[(f64, f64)],
    summary: &SpectraSummary,
) -> anyhow::Result<()> {
    let views = build_views(hist, summary);
    let mut template = include_str!("../../assets/spectra.tex").to_string();
    super::common::replace_section(&mut template, "%SPECTRA_BEGIN", "%SPECTRA_END", &views)?;
    w.write_all(template.as_bytes())?;
    Ok(())
}

/// Generate the standard spectra view: observed histogram, full model,
/// error region, and k-mer peak markers (linear coverage vs frequency).
fn build_views(hist: &[(f64, f64)], summary: &SpectraSummary) -> String {
    let params = &summary.params;
    let kcovfloor = ((params.kmercov - 2.0 * summary.se_kmercov).floor().max(1.0)) as usize;
    // R `error_xcutoff_ind = tail(which(x <= error_xcutoff), 1)`.
    let err_cut = hist
        .iter()
        .rposition(|&(x, _)| x <= kcovfloor as f64)
        .map_or(1, |i| i + 1);

    // Error residuals (count space), truncated at the first value < 1.
    let mut error_kmers = vec![0.0f64; hist.len()];
    let mut first_zero = err_cut;
    for i in 0..err_cut {
        let (x, count) = hist[i];
        let pred = params.length * predict(summary.p, 1, params, summary.k, x);
        let e = count - pred;
        if e < 1.0 {
            first_zero = i;
            break;
        }
        error_kmers[i] = e;
    }

    // Observed and model series.
    let mut linear_obs = String::new();
    let mut linear_model = String::new();
    let mut errors = String::new();
    for (i, &(x, count)) in hist.iter().enumerate() {
        let model = params.length * predict(summary.p, 1, params, summary.k, x);
        linear_obs.push_str(&format!("({:.1},{:.0}) ", x, count));
        linear_model.push_str(&format!("({:.1},{:.4}) ", x, model));
        if i < first_zero {
            errors.push_str(&format!("({:.1},{:.4}) ", x, error_kmers[i]));
        }
    }

    let title = format!(
        "len {:.0} bp, uniq {:.2}\\%, kcov {:.2}, err {:.2}\\%, dup {:.3} \\\\ k {}, p {}",
        summary.len,
        summary.unique * 100.0,
        params.kmercov,
        summary.error_rate * 100.0,
        params.d,
        summary.k,
        summary.p
    );

    // Peak markers at kmercov * 1..2p.
    let mut peaks = String::new();
    let ymax = hist.iter().map(|&(_, c)| c).fold(0.0f64, f64::max).max(1.0);
    for i in 1..=2 * summary.p {
        let x = params.kmercov * i as f64;
        if x <= ymax * 1.5 {
            peaks.push_str(&format!("({:.2},0) ({:.2},{:.0}) ", x, x, ymax));
        }
    }

    format!(
        r###"\begin{{tikzpicture}}
  \begin{{axis}}[
    title={{{title}}},
    title style={{align=center}},
    xlabel={{Coverage}}, ylabel={{Frequency}},
    axis background/.style={{fill=lightgray!20}},
    width=9cm, height=7cm,
    legend pos=north east,
  ]
    \addplot[color=blue!70!black, ycomb, mark=none, line width=0.8pt] coordinates {{ {linear_obs} }};
    \addplot[black, very thick] coordinates {{ {linear_model} }};
    \addplot[orange!80!black, very thick] coordinates {{ {errors} }};
    \addplot[dashed, black] coordinates {{ {peaks} }};
    \legend{{observed, full model, errors, kmer peaks}}
  \end{{axis}}
\end{{tikzpicture}}

"###,
        title = title,
        linear_obs = linear_obs,
        linear_model = linear_model,
        errors = errors,
        peaks = peaks,
    )
}
