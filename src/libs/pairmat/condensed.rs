/// A condensed distance matrix (upper triangle only, no diagonal).
///
/// Stores only the upper triangular part of a symmetric matrix, reducing memory usage
/// from N^2 to N(N-1)/2. This format is required by some hierarchical clustering algorithms
/// (like `kodama`'s linkage) and is the core storage engine for `NamedMatrix`.
///
/// # Storage Layout
/// For N=3, indices (0,1), (0,2), (1,2) are stored at 0, 1, 2.
/// The diagonal is implicitly 0.0. The matrix is symmetric.
///
/// # Examples
/// ```
/// # use pgr::libs::pairmat::CondensedMatrix;
/// let mut m = CondensedMatrix::new(3);
/// m.set(0, 1, 0.5);
/// m.set(0, 2, 0.8);
/// m.set(1, 2, 0.3);
///
/// assert_eq!(m.get(0, 1), 0.5);
/// assert_eq!(m.get(1, 0), 0.5); // Symmetric
/// assert_eq!(m.get(0, 0), 0.0); // Diagonal is always 0
/// ```
#[derive(Debug, Clone)]
pub struct CondensedMatrix {
    size: usize,
    data: Vec<f32>,
}

impl CondensedMatrix {
    /// Create a new condensed matrix of size N x N.
    pub fn new(size: usize) -> Self {
        let len = if size == 0 { 0 } else { size * (size - 1) / 2 };
        Self {
            size,
            data: vec![0.0; len],
        }
    }

    /// Create from existing data vector.
    ///
    /// # Panics
    /// Panics if data length doesn't match size*(size-1)/2.
    pub fn from_vec(size: usize, data: Vec<f32>) -> Self {
        let expected = if size == 0 { 0 } else { size * (size - 1) / 2 };
        assert_eq!(
            data.len(),
            expected,
            "Data length {} does not match expected length {} for size {}",
            data.len(),
            expected,
            size
        );
        Self { size, data }
    }

    pub fn size(&self) -> usize {
        self.size
    }

    /// Get the underlying data vector.
    pub fn data(&self) -> &[f32] {
        &self.data
    }

    /// Get the underlying data vector as mutable slice.
    pub fn data_mut(&mut self) -> &mut [f32] {
        &mut self.data
    }

    /// Get value at (row, col).
    /// Returns 0.0 if row == col.
    pub fn get(&self, row: usize, col: usize) -> f32 {
        if row == col {
            0.0
        } else if row < col {
            self.data[get_condensed_index(self.size, row, col)]
        } else {
            self.data[get_condensed_index(self.size, col, row)]
        }
    }

    /// Set value at (row, col).
    /// Does nothing if row == col.
    pub fn set(&mut self, row: usize, col: usize, value: f32) {
        if row != col {
            let idx = if row < col {
                get_condensed_index(self.size, row, col)
            } else {
                get_condensed_index(self.size, col, row)
            };
            self.data[idx] = value;
        }
    }
}

/// Convert row, col (where row < col) to linear index in condensed array (upper triangle).
///
/// Based on formula: k = N*row - row*(row+1)/2 + col - row - 1
#[inline]
pub fn get_condensed_index(size: usize, row: usize, col: usize) -> usize {
    debug_assert!(row < col);
    debug_assert!(col < size);
    size * row - (row * (row + 1)) / 2 + col - row - 1
}
