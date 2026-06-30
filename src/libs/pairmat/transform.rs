use super::named::NamedMatrix;

/// Apply mathematical transformations to a matrix element-wise.
///
/// Supports: linear, inv-linear, log, exp, square, sqrt.
/// When `normalize` is true, off-diagonal values are divided by `sqrt(d_i * d_j)`
/// and diagonal values are normalized to 1.0 (or 0.0 if the original diag <= 1e-9).
pub fn transform_matrix(
    matrix: &NamedMatrix,
    method: &str,
    max_val: f32,
    scale: f32,
    offset: f32,
    normalize: bool,
) -> anyhow::Result<NamedMatrix> {
    let mut result = matrix.clone();
    let size = result.size();

    // Get original diagonals (used for normalize and for transforming diagonal elements)
    let diags: Vec<f32> = result.get_diags().cloned().unwrap_or_default();
    let has_diags = !diags.is_empty();

    // Warn if normalize is requested but diagonals are missing or all zero
    if normalize {
        if diags.is_empty() {
            log::warn!(
                "--normalize requested but no diagonal values found. Result will be Inf/NaN."
            );
        } else {
            let max_diag = diags.iter().fold(0.0f32, |a, &b| a.max(b));
            if max_diag == 0.0 {
                log::warn!("--normalize requested but all diagonal values are 0.0. Result will be Inf/NaN.");
            }
        }
    }

    // Transform off-diagonal elements (upper triangle)
    for i in 0..size {
        for j in (i + 1)..size {
            let mut val = result.get(i, j);

            // 1. Normalize
            if normalize {
                let d_i = diags[i];
                let d_j = diags[j];
                if d_i > 1e-9 && d_j > 1e-9 {
                    val /= (d_i * d_j).sqrt();
                } else {
                    val = 0.0;
                }
            }

            // 2. Transform
            val = match method {
                "linear" => val * scale + offset,
                "inv-linear" => max_val - val,
                "log" => {
                    if val > 0.0 {
                        -val.ln()
                    } else {
                        1000.0
                    }
                }
                "exp" => (-val).exp(),
                "square" => val * val,
                "sqrt" => {
                    if val >= 0.0 {
                        val.sqrt()
                    } else {
                        0.0
                    }
                }
                _ => val,
            };

            result.set(i, j, val);
        }
    }

    // Transform diagonal elements.
    // Normalize sets d to 1.0 (if original d > 1e-9) or 0.0, matching off-diagonal behavior
    // where x_norm(i,i) = x(i,i) / sqrt(x(i,i)*x(i,i)) = 1.0.
    let mut new_diags = vec![0.0; size];
    for i in 0..size {
        let mut d = if has_diags { diags[i] } else { 0.0 };
        if normalize {
            d = if d > 1e-9 { 1.0 } else { 0.0 };
        }
        d = match method {
            "linear" => d * scale + offset,
            "inv-linear" => max_val - d,
            "log" => {
                if d > 0.0 {
                    -d.ln()
                } else {
                    0.0
                }
            }
            "exp" => (-d).exp(),
            "square" => d * d,
            "sqrt" => {
                if d >= 0.0 {
                    d.sqrt()
                } else {
                    0.0
                }
            }
            _ => d,
        };
        new_diags[i] = d;
    }
    result.set_diags(new_diags);

    Ok(result)
}
