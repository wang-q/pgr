use std::cmp;

/// A gap cost calculator using linear interpolation for efficient scoring.
///
/// It uses pre-calculated tables for small gap sizes and interpolation for larger ones.
/// Separate costs are maintained for query gaps, target gaps, and simultaneous gaps (both).
#[derive(Clone, Debug)]
pub struct GapCalc {
    small_size: usize,
    q_small: Vec<i32>,
    t_small: Vec<i32>,
    b_small: Vec<i32>,
    long_pos: Vec<i32>,
    q_long: Vec<f64>,
    t_long: Vec<f64>,
    b_long: Vec<f64>,
    // Last params for extrapolation
    q_last_pos: i32,
    q_last_pos_val: f64,
    q_last_slope: f64,
    t_last_pos: i32,
    t_last_pos_val: f64,
    t_last_slope: f64,
    b_last_pos: i32,
    b_last_pos_val: f64,
    b_last_slope: f64,
}

impl GapCalc {
    /// Creates a standard "medium" gap calculator (suitable for mouse/human).
    pub fn medium() -> Self {
        // "medium" (mouse/human)
        // position: 1, 2, 3, 11, 111, 2111, 12111, 32111, 72111, 152111, 252111
        let pos = vec![1, 2, 3, 11, 111, 2111, 12111, 32111, 72111, 152111, 252111];
        let q_gap = vec![
            325.0, 360.0, 400.0, 450.0, 600.0, 1100.0, 3600.0, 7600.0, 15600.0, 31600.0, 56600.0,
        ];
        // tGap same as qGap
        let t_gap = q_gap.clone();
        let b_gap = vec![
            625.0, 660.0, 700.0, 750.0, 900.0, 1400.0, 4000.0, 8000.0, 16000.0, 32000.0, 57000.0,
        ];

        Self::new(pos, q_gap, t_gap, b_gap)
    }

    /// Creates a "loose" gap calculator (suitable for distant species like chicken/human).
    pub fn loose() -> Self {
        // "loose" (chicken/human)
        // position: 1, 2, 3, 11, 111, 2111, 12111, 32111, 72111, 152111, 252111
        let pos = vec![1, 2, 3, 11, 111, 2111, 12111, 32111, 72111, 152111, 252111];
        let q_gap = vec![
            350.0, 425.0, 450.0, 600.0, 900.0, 2900.0, 22900.0, 57900.0, 117900.0, 217900.0,
            317900.0,
        ];
        let b_gap = vec![
            750.0, 825.0, 850.0, 1000.0, 1300.0, 3300.0, 23300.0, 58300.0, 118300.0, 218300.0,
            318300.0,
        ];
        let t_gap = q_gap.clone();

        Self::new(pos, q_gap, t_gap, b_gap)
    }

    /// Creates a new `GapCalc` with affine gap costs.
    ///
    /// # Arguments
    ///
    /// * `open` - Gap open cost.
    /// * `extend` - Gap extension cost.
    pub fn affine(open: i32, extend: i32) -> Self {
        let pos = vec![1, 2, 3, 11, 111, 2111, 12111, 32111, 72111, 152111, 252111];

        let calc_cost = |len: i32| -> f64 {
            if len <= 0 {
                0.0
            } else {
                (open + extend * len) as f64
            }
        };

        let q_gap: Vec<f64> = pos.iter().map(|&x| calc_cost(x)).collect();
        let t_gap = q_gap.clone();

        // For simultaneous gaps, we can use max(dq, dt) logic in calc(),
        // but here we just need a table.
        // In affine mode, calc() handles sim gaps by taking max(dq, dt) and looking up single gap cost.
        // So b_gap table should be same as q_gap/t_gap.
        let b_gap = q_gap.clone();

        Self::new(pos, q_gap, t_gap, b_gap)
    }

    /// Creates a new `GapCalc` with custom cost tables.
    ///
    /// # Arguments
    ///
    /// * `pos` - Positions (gap sizes) for which costs are defined.
    /// * `q_vals` - Costs for gaps in query sequence.
    /// * `t_vals` - Costs for gaps in target sequence.
    /// * `b_vals` - Costs for gaps in both sequences (simultaneous).
    pub fn new(pos: Vec<i32>, q_vals: Vec<f64>, t_vals: Vec<f64>, b_vals: Vec<f64>) -> Self {
        let small_size = 111;
        let mut q_small = vec![0; small_size];
        let mut t_small = vec![0; small_size];
        let mut b_small = vec![0; small_size];

        for i in 1..small_size {
            q_small[i] = Self::interpolate(i as i32, &pos, &q_vals) as i32;
            t_small[i] = Self::interpolate(i as i32, &pos, &t_vals) as i32;
            b_small[i] = Self::interpolate(i as i32, &pos, &b_vals) as i32;
        }

        let start_long = pos
            .iter()
            .position(|&x| x == small_size as i32)
            .unwrap_or(0);

        let long_pos = pos[start_long..].to_vec();
        let q_long = q_vals[start_long..].to_vec();
        let t_long = t_vals[start_long..].to_vec();
        let b_long = b_vals[start_long..].to_vec();

        let n = long_pos.len();
        let q_last_pos = long_pos[n - 1];
        let q_last_pos_val = q_long[n - 1];
        let q_last_slope =
            (q_long[n - 1] - q_long[n - 2]) / (long_pos[n - 1] - long_pos[n - 2]) as f64;

        let t_last_pos = long_pos[n - 1];
        let t_last_pos_val = t_long[n - 1];
        let t_last_slope =
            (t_long[n - 1] - t_long[n - 2]) / (long_pos[n - 1] - long_pos[n - 2]) as f64;

        let b_last_pos = long_pos[n - 1];
        let b_last_pos_val = b_long[n - 1];
        let b_last_slope =
            (b_long[n - 1] - b_long[n - 2]) / (long_pos[n - 1] - long_pos[n - 2]) as f64;

        GapCalc {
            small_size,
            q_small,
            t_small,
            b_small,
            long_pos,
            q_long,
            t_long,
            b_long,
            q_last_pos,
            q_last_pos_val,
            q_last_slope,
            t_last_pos,
            t_last_pos_val,
            t_last_slope,
            b_last_pos,
            b_last_pos_val,
            b_last_slope,
        }
    }

    fn interpolate(x: i32, s: &[i32], v: &[f64]) -> f64 {
        for i in 0..s.len() {
            if x == s[i] {
                return v[i];
            } else if x < s[i] {
                if i == 0 {
                    return v[0];
                }
                let ds = s[i] - s[i - 1];
                let dv = v[i] - v[i - 1];
                return v[i - 1] + dv * (x - s[i - 1]) as f64 / ds as f64;
            }
        }
        let n = s.len();
        let ds = s[n - 1] - s[n - 2];
        let dv = v[n - 1] - v[n - 2];
        v[n - 2] + dv * (x - s[n - 2]) as f64 / ds as f64
    }

    /// Calculates the gap cost for a given distance in query (`dq`) and target (`dt`).
    pub fn calc(&self, dq: i32, dt: i32) -> i32 {
        let dt = if dt < 0 { 0 } else { dt };
        let dq = if dq < 0 { 0 } else { dq };

        if dt == 0 {
            if (dq as usize) < self.small_size {
                self.q_small[dq as usize]
            } else if dq >= self.q_last_pos {
                let cost = self.q_last_pos_val + self.q_last_slope * (dq - self.q_last_pos) as f64;
                cost as i32
            } else {
                Self::interpolate(dq, &self.long_pos, &self.q_long) as i32
            }
        } else if dq == 0 {
            if (dt as usize) < self.small_size {
                self.t_small[dt as usize]
            } else if dt >= self.t_last_pos {
                let cost = self.t_last_pos_val + self.t_last_slope * (dt - self.t_last_pos) as f64;
                cost as i32
            } else {
                Self::interpolate(dt, &self.long_pos, &self.t_long) as i32
            }
        } else {
            // For simultaneous gaps, we use max(dq, dt) to determine the cost
            let both = cmp::max(dq, dt);
            if (both as usize) < self.small_size {
                self.b_small[both as usize]
            } else if both >= self.b_last_pos {
                let cost =
                    self.b_last_pos_val + self.b_last_slope * (both - self.b_last_pos) as f64;
                cost as i32
            } else {
                Self::interpolate(both, &self.long_pos, &self.b_long) as i32
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gap_calc_medium() {
        let calc = GapCalc::medium();

        // Test small values (should be in small table)
        // pos: 1 -> 325.0
        assert_eq!(calc.calc(1, 0), 325);
        assert_eq!(calc.calc(0, 1), 325);
    }
}
