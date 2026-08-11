//! BBMerge overlap neural network (bbmerge.bbnet) inference.

use anyhow::{Context, Result};
use std::path::Path;

/// Activation functions used by the BBTools `CellNet` format.
#[derive(Debug, Clone, Copy)]
enum Activation {
    Sigmoid,
    Tanh,
    MSig,
    RSLog,
    Linear,
}

/// One dense-layer neuron.
#[derive(Debug, Clone)]
struct Cell {
    function: Activation,
    bias: f32,
    weights: Vec<f32>,
}

/// A parsed dense `CellNet` (BBMerge overlap filter).
#[derive(Debug, Clone)]
pub struct CellNet {
    layers: Vec<Vec<Cell>>,
    /// Classification threshold (`##ctf`); scores below it are rejected.
    pub cutoff: f32,
}

impl CellNet {
    /// Parses a `.bbnet` file (dense, concise format).
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let text = std::fs::read_to_string(path.as_ref())
            .with_context(|| format!("failed to read {}", path.as_ref().display()))?;
        Self::parse(&text)
    }

    /// Parses `.bbnet` text (dense, concise format).
    pub fn parse(text: &str) -> Result<Self> {
        let mut dims: Vec<usize> = Vec::new();
        let mut cutoff = 0.5f32;
        let mut cells: Vec<(usize, Activation, f32, Vec<f32>)> = Vec::new();

        for line in text.lines() {
            let line = line.trim();
            if line.starts_with("#dims") {
                dims = line
                    .split_whitespace()
                    .skip(1)
                    .map(|s| s.parse::<usize>())
                    .collect::<std::result::Result<Vec<_>, _>>()
                    .context("bad #dims")?;
            } else if let Some(rest) = line.strip_prefix("##ctf") {
                cutoff = rest.trim().parse::<f32>().context("bad ##ctf")?;
            } else if line.starts_with('C') || line.starts_with('W') {
                let mut parts = line.split_whitespace();
                let id = parts
                    .next()
                    .unwrap()
                    .trim_start_matches(['C', 'W'])
                    .parse::<usize>()
                    .context("bad cell id")?;
                let fname = parts.next().context("missing function")?;
                let function = match fname {
                    "SIG" | "SIGMOID" => Activation::Sigmoid,
                    "TANH" | "HYPERBOLICTANGENT" => Activation::Tanh,
                    "MSIG" | "MIRROREDSIGMOID" => Activation::MSig,
                    "RSLOG" | "ROTATIONALLYSYMMETRICLOGARITHM" => Activation::RSLog,
                    "LINEAR" => Activation::Linear,
                    _ => anyhow::bail!("unknown net activation {fname}"),
                };
                let bias = parts.next().context("missing bias")?.parse::<f32>()?;
                let weights = parts
                    .map(|s| s.parse::<f32>())
                    .collect::<std::result::Result<Vec<_>, _>>()
                    .context("bad weights")?;
                cells.push((id, function, bias, weights));
            }
        }
        anyhow::ensure!(dims.len() >= 2, "missing or bad #dims");

        // Cells are numbered across layers starting at 1 (index 0 reserved):
        // layer l holds dims[l] cells with ids [prefix, prefix+dims[l]).
        let mut layers: Vec<Vec<Cell>> = Vec::with_capacity(dims.len());
        let mut start = 1usize;
        for (l, &size) in dims.iter().enumerate() {
            if l == 0 {
                // Input layer: no weights.
                layers.push(Vec::new());
                start += size; // 1 + dims[0]
                continue;
            }
            let mut layer = Vec::with_capacity(size);
            for id in start..start + size {
                let (fid, function, bias, weights) = cells
                    .iter()
                    .find(|(cid, _, _, _)| *cid == id)
                    .with_context(|| format!("missing cell {id}"))?
                    .clone();
                let _ = fid;
                layer.push(Cell {
                    function,
                    bias,
                    weights,
                });
            }
            layers.push(layer);
            start += size;
        }
        Ok(CellNet { layers, cutoff })
    }

    /// Runs the 23-feature vector through the network.
    pub fn feed_forward(&self, input: &[f32]) -> f32 {
        let mut values = vec![input.to_vec()];
        for layer in self.layers.iter().skip(1) {
            let prev = values.last().unwrap();
            let mut out = Vec::with_capacity(layer.len());
            for c in layer {
                // BBTools Vector.fma: 256-bit SIMD (8 lanes of fused
                // multiply-add) when the previous layer has >= 32 cells,
                // otherwise a plain scalar sum. The lane reduction order
                // mirrors the JDK FloatVector reduceLanes for 256-bit.
                let sum = c.bias
                    + if prev.len() >= 32 {
                        simd_fma_256(prev, &c.weights)
                    } else {
                        scalar_dot(prev, &c.weights)
                    };
                out.push(activate(c.function, sum));
            }
            values.push(out);
        }
        *values.last().unwrap().first().unwrap()
    }
}

/// Plain f32 dot product (BBTools scalar fallback).
fn scalar_dot(a: &[f32], b: &[f32]) -> f32 {
    let mut c = 0f32;
    for i in 0..a.len() {
        c += a[i] * b[i];
    }
    c
}

/// 256-bit SIMD dot product with fused multiply-add lanes.
fn simd_fma_256(a: &[f32], b: &[f32]) -> f32 {
    let limit = a.len() - a.len() % 8;
    let mut lanes = [0f32; 8];
    let mut i = 0;
    while i < limit {
        for l in 0..8 {
            lanes[l] = a[i + l].mul_add(b[i + l], lanes[l]);
        }
        i += 8;
    }
    // JDK FloatVector reduceLanes(ADD): balanced binary tree over the 8 lanes.
    let mut c = 0f32;
    for l in 0..4 {
        c += lanes[l] + lanes[l + 4];
    }
    for j in i..a.len() {
        c += a[j] * b[j];
    }
    c
}

/// Applies a BBTools activation (double precision, cast to f32).
fn activate(f: Activation, x: f32) -> f32 {
    let x = x as f64;
    let y = match f {
        Activation::Sigmoid => 1.0 / (1.0 + (-x).exp()),
        Activation::Tanh => {
            if x < -20.0 {
                -1.0
            } else if x > 20.0 {
                1.0
            } else {
                let ex = x.exp();
                let emx = (-x).exp();
                (ex - emx) / (ex + emx)
            }
        }
        Activation::MSig => {
            let offset = 5.0f64;
            let xmult = 2.0f64;
            // MSIG_Y_MULT = 1/sigmoid(offset) = 1+exp(-offset).
            let ymult = 1.0 + (-offset).exp();
            let y = if x < 0.0 {
                1.0 / (1.0 + (-(xmult * x + offset)).exp())
            } else {
                1.0 / (1.0 + (xmult * x - offset).exp())
            };
            ymult * y
        }
        Activation::RSLog => {
            if x < 0.0 {
                -(-x + 1.0).ln()
            } else {
                (x + 1.0).ln()
            }
        }
        Activation::Linear => x,
    };
    y as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_and_forward_minimal() {
        // 2 inputs -> 2 hidden -> 1 output.
        let text = "\
##bbnet
#dims 2 2 1
##ctf 0.5
C3 TANH 0.0 1.0 1.0
C4 SIG 0.0 -1.0 -1.0
C5 SIG 0.0 1.0 1.0
";
        let net = CellNet::parse(text).unwrap();
        let out = net.feed_forward(&[1.0, 1.0]);
        assert!(out > 0.5);
    }
}
