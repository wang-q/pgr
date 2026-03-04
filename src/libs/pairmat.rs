//! A *symmetric* scoring matrix to be used for clustering.
use std::collections::HashMap;
use std::io::BufRead;

#[derive(Debug, Clone)]
pub struct ScoringMatrix<T> {
    size: Option<usize>,
    same: Option<T>,
    missing: Option<T>,
    data: HashMap<(usize, usize), T>,
}

impl<T> ScoringMatrix<T>
where
    T: Default + Copy,
{
    /// Creates a new empty matrix with default values.
    ///
    /// ```
    /// # use pgr::libs::pairmat::ScoringMatrix;
    /// let matrix: ScoringMatrix<i32> = ScoringMatrix::new();
    /// assert_eq!(matrix.get(0, 0), 0);  // Using T::default()
    /// ```
    pub fn new() -> Self {
        ScoringMatrix {
            size: None,
            same: None,
            missing: None,
            data: HashMap::new(),
        }
    }

    /// Creates a new matrix with specified default values.
    ///
    /// ```
    /// # use pgr::libs::pairmat::ScoringMatrix;
    /// let matrix = ScoringMatrix::with_defaults(0.0, -1.0);
    /// assert_eq!(matrix.get(0, 0), 0.0);    // same value
    /// assert_eq!(matrix.get(0, 1), -1.0);   // missing value
    /// ```
    pub fn with_defaults(same: T, missing: T) -> Self {
        ScoringMatrix {
            size: None,
            same: Some(same),
            missing: Some(missing),
            data: HashMap::new(),
        }
    }

    /// Creates a new matrix with specified size and default values.
    ///
    /// ```
    /// # use pgr::libs::pairmat::ScoringMatrix;
    /// let matrix = ScoringMatrix::with_size_and_defaults(3, 1.0, 0.0);
    /// assert_eq!(matrix.size(), 3);
    /// assert_eq!(matrix.get(0, 0), 1.0);    // same value
    /// assert_eq!(matrix.get(0, 1), 0.0);    // missing value
    /// ```
    pub fn with_size_and_defaults(size: usize, same: T, missing: T) -> Self {
        ScoringMatrix {
            size: Some(size),
            same: Some(same),
            missing: Some(missing),
            data: HashMap::new(),
        }
    }

    pub fn with_size(size: usize) -> Self {
        ScoringMatrix {
            size: Some(size),
            same: None,
            missing: None,
            data: HashMap::new(),
        }
    }

    pub fn size(&self) -> usize {
        self.size.unwrap_or_else(|| {
            self.data
                .keys()
                .map(|&(i, j)| i.max(j) + 1)
                .max()
                .unwrap_or(0)
        })
    }

    /// Sets a fixed size for the matrix
    pub fn set_size(&mut self, size: usize) {
        self.size = Some(size);
    }

    /// Returns the value of the given cell.
    ///
    /// ```
    /// # use pgr::libs::pairmat::ScoringMatrix;
    /// let mut m = ScoringMatrix::with_size_and_defaults(5, 0, 1);
    /// m.set(1, 2, 42);
    /// assert_eq!(m.get(1, 2), 42);
    /// assert_eq!(m.get(2, 1), 42);
    /// assert_eq!(m.get(3, 3), 0);
    /// assert_eq!(m.get(1, 3), 1);
    /// ```
    pub fn set(&mut self, row: usize, col: usize, value: T) {
        if row <= col {
            self.data.insert((row, col), value);
        } else {
            self.data.insert((col, row), value);
        }
    }

    pub fn get(&self, row: usize, col: usize) -> T {
        if row == col {
            self.data
                .get(&(row, col))
                .copied()
                .unwrap_or_else(|| self.same.unwrap_or(T::default()))
        } else {
            let (r, c) = if row < col { (row, col) } else { (col, row) };
            self.data
                .get(&(r, c))
                .copied()
                .unwrap_or_else(|| self.missing.unwrap_or(T::default()))
        }
    }
}

// Add a separate implementation for f32 specifically for from_pair_scores
impl ScoringMatrix<f32> {
    pub fn from_pair_scores(infile: &str, same: f32, missing: f32) -> (Self, Vec<String>) {
        let mut names = indexmap::IndexSet::new();
        let mut matrix = Self::with_defaults(same, missing);

        let reader = crate::reader(infile);
        for line in reader.lines().map_while(Result::ok) {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 3 {
                let n1 = fields[0].to_string();
                let n2 = fields[1].to_string();
                let score = fields[2].parse().unwrap();

                names.insert(n1.clone());
                names.insert(n2.clone());

                matrix.set(
                    names.get_index_of(&n1).unwrap(),
                    names.get_index_of(&n2).unwrap(),
                    score,
                );
            }
        }

        matrix.set_size(names.len());
        (matrix, names.into_iter().collect())
    }
}

/// A named matrix for storing pairwise distances/scores with sequence names
///
/// # Examples
///
/// ```
/// # use pgr::libs::pairmat::NamedMatrix;
/// let names = vec!["seq1".to_string(), "seq2".to_string(), "seq3".to_string()];
/// let mut matrix = NamedMatrix::new(names);
///
/// // Set some values
/// matrix.set(0, 1, 0.5);
/// matrix.set(0, 2, 0.7);
/// matrix.set(1, 2, 0.3);
///
/// // Get values
/// assert_eq!(matrix.size(), 3);
/// assert_eq!(matrix.get(0, 1), 0.5);
/// assert_eq!(matrix.get(1, 0), 0.5);  // Symmetric matrix
/// ```
#[derive(Debug, Clone)]
pub struct NamedMatrix {
    size: usize,
    names: indexmap::IndexMap<String, usize>,
    values: Vec<f32>,
}

impl NamedMatrix {
    pub fn new(names: Vec<String>) -> Self {
        let size = names.len();
        let values = vec![f32::default(); size * size];
        let names: indexmap::IndexMap<_, _> = names
            .into_iter()
            .enumerate()
            .map(|(i, name)| (name, i))
            .collect();

        NamedMatrix {
            size,
            names,
            values,
        }
    }

    pub fn size(&self) -> usize {
        self.size
    }

    pub fn get(&self, row: usize, col: usize) -> f32 {
        self.values[row * self.size + col]
    }

    pub fn set(&mut self, row: usize, col: usize, value: f32) {
        self.values[row * self.size + col] = value;
        // Set the symmetric position as it's a symmetric matrix
        if row != col {
            self.values[col * self.size + row] = value;
        }
    }

    pub fn get_names(&self) -> Vec<&String> {
        self.names.keys().collect()
    }

    /// Get matrix value by sequence names
    ///
    /// ```
    /// # use pgr::libs::pairmat::NamedMatrix;
    /// let names = vec!["seq1".to_string(), "seq2".to_string()];
    /// let mut matrix = NamedMatrix::new(names);
    /// matrix.set(0, 1, 0.5);
    ///
    /// assert_eq!(matrix.get_by_name("seq1", "seq2"), Some(0.5));
    /// assert_eq!(matrix.get_by_name("seq1", "seq3"), None);  // Non-existent name
    /// ```
    pub fn get_by_name(&self, name1: &str, name2: &str) -> Option<f32> {
        let i = self.names.get(name1)?;
        let j = self.names.get(name2)?;
        Some(self.get(*i, *j))
    }

    /// Set matrix value by sequence names
    ///
    /// ```
    /// # use pgr::libs::pairmat::NamedMatrix;
    /// let names = vec!["seq1".to_string(), "seq2".to_string()];
    /// let mut matrix = NamedMatrix::new(names);
    ///
    /// assert!(matrix.set_by_name("seq1", "seq2", 0.5).is_ok());
    /// assert_eq!(matrix.get_by_name("seq1", "seq2"), Some(0.5));
    /// assert!(matrix.set_by_name("seq1", "seq3", 0.5).is_err());  // Non-existent name
    /// ```
    pub fn set_by_name(&mut self, name1: &str, name2: &str, value: f32) -> Result<(), String> {
        match (self.names.get(name1), self.names.get(name2)) {
            (Some(&i), Some(&j)) => {
                self.set(i, j, value);
                Ok(())
            }
            (None, _) => Err(format!("Name not found: {}", name1)),
            (_, None) => Err(format!("Name not found: {}", name2)),
        }
    }

    pub fn from_pair_scores(infile: &str, same: f32, missing: f32) -> Self {
        let (scoring_matrix, index_name) = ScoringMatrix::from_pair_scores(infile, same, missing);
        let size = index_name.len();

        // Create NamedMatrix from ScoringMatrix
        let mut values = vec![f32::default(); size * size];
        for i in 0..size {
            for j in 0..size {
                values[i * size + j] = scoring_matrix.get(i, j);
            }
        }

        // Convert index_name to IndexMap
        let names: indexmap::IndexMap<_, _> = index_name
            .into_iter()
            .enumerate()
            .map(|(i, name)| (name, i))
            .collect();

        Self {
            size,
            names,
            values,
        }
    }

    /// Creates a new matrix from a relaxed PHYLIP format file
    ///
    /// ```no_run
    /// # use pgr::libs::pairmat::NamedMatrix;
    /// let matrix = NamedMatrix::from_relaxed_phylip("input.phy");
    /// ```
    pub fn from_relaxed_phylip(infile: &str) -> Self {
        let mut names = Vec::new();
        let mut values = Vec::new();

        let reader = crate::reader(infile);
        let mut lines = reader.lines();

        // Skip the optional sequence count line
        if let Some(Ok(line)) = lines.next() {
            if line.trim().parse::<usize>().is_err() {
                // If first line is not a number, treat it as a data line
                Self::process_phylip_line(&line, &mut names, &mut values);
            }
        }

        // Process remaining lines
        for line in lines.map_while(Result::ok) {
            Self::process_phylip_line(&line, &mut names, &mut values);
        }

        let size = names.len();
        let mut matrix = Self::new(names);

        // Fill the matrix
        for i in 0..size {
            for j in 0..=i {
                let value = values[i * (i + 1) / 2 + j];
                matrix.set(i, j, value);
            }
        }

        matrix
    }

    fn process_phylip_line(line: &str, names: &mut Vec<String>, values: &mut Vec<f32>) {
        let parts: Vec<&str> = line.trim().split_whitespace().collect();
        if !parts.is_empty() {
            let name = parts[0].to_string();
            names.push(name);

            // Read lower-triangle distances
            let distances: Vec<f32> = parts[1..=names.len()]
                .iter()
                .map(|&s| s.parse().unwrap())
                .collect();

            values.extend(distances);
        }
    }
}

/// A condensed distance matrix (upper triangle only, no diagonal).
///
/// Stores only the upper triangular part of a symmetric matrix, reducing memory usage
/// from N^2 to N(N-1)/2. This format is required by some hierarchical clustering algorithms
/// (like `kodama`'s linkage).
///
/// # Storage Layout
/// For N=3, indices (0,1), (0,2), (1,2) are stored at 0, 1, 2.
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
    data: Vec<f32>, // length = N*(N-1)/2
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

    /// Convert row, col (where row < col) to linear index in condensed array.
    /// Based on formula: k = N*row - row*(row+1)/2 + col - row - 1
    #[inline]
    fn index(&self, row: usize, col: usize) -> usize {
        debug_assert!(row < col);
        debug_assert!(col < self.size);
        // Standard formula for row-major upper triangle without diagonal
        // sum(N-1-i for i in 0..row) + (col - row - 1)
        // = (N-1)*row - row*(row-1)/2 + col - row - 1
        // = N*row - row - row^2/2 + row/2 + col - row - 1
        // = N*row - row*(row+1)/2 + col - row - 1 (simplifies to this)
        
        // Let's use the iterative sum logic which is safer to reason about:
        // offset for row `r` is sum_{i=0}^{r-1} (N - 1 - i)
        // number of elements in row `i` is (N - 1 - i)
        
        // Using the closed form:
        // offset = r*N - r*(r+1)/2
        // index = offset + (c - r - 1)
        
        self.size * row - (row * (row + 1)) / 2 + col - row - 1
    }

    /// Get value at (row, col).
    /// Returns 0.0 if row == col.
    pub fn get(&self, row: usize, col: usize) -> f32 {
        if row == col {
            0.0
        } else if row < col {
            self.data[self.index(row, col)]
        } else {
            self.data[self.index(col, row)]
        }
    }

    /// Set value at (row, col).
    /// Does nothing if row == col.
    pub fn set(&mut self, row: usize, col: usize, value: f32) {
        if row != col {
            let idx = if row < col {
                self.index(row, col)
            } else {
                self.index(col, row)
            };
            self.data[idx] = value;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_condensed_matrix_indexing() {
        // N=4
        // (0,1) -> 0
        // (0,2) -> 1
        // (0,3) -> 2
        // (1,2) -> 3
        // (1,3) -> 4
        // (2,3) -> 5
        let m = CondensedMatrix::new(4);
        assert_eq!(m.index(0, 1), 0);
        assert_eq!(m.index(0, 2), 1);
        assert_eq!(m.index(0, 3), 2);
        assert_eq!(m.index(1, 2), 3);
        assert_eq!(m.index(1, 3), 4);
        assert_eq!(m.index(2, 3), 5);
    }

    #[test]
    fn test_condensed_matrix_rw() {
        let mut m = CondensedMatrix::new(3);
        m.set(0, 1, 1.0);
        m.set(2, 0, 2.0); // set (0,2) via swap
        m.set(1, 2, 3.0);

        assert_eq!(m.get(0, 1), 1.0);
        assert_eq!(m.get(1, 0), 1.0);
        assert_eq!(m.get(0, 2), 2.0);
        assert_eq!(m.get(2, 0), 2.0);
        assert_eq!(m.get(1, 2), 3.0);
        assert_eq!(m.get(0, 0), 0.0);
    }

    #[test]
    fn test_scoring_matrix_basic() {
        let mut m = ScoringMatrix::with_defaults(0.0, -1.0);
        m.set(0, 1, 5.0);
        m.set(2, 1, 10.0);

        // Check set values (symmetric)
        assert_eq!(m.get(0, 1), 5.0);
        assert_eq!(m.get(1, 0), 5.0);
        assert_eq!(m.get(1, 2), 10.0);
        
        // Check diagonal default
        assert_eq!(m.get(0, 0), 0.0);
        assert_eq!(m.get(3, 3), 0.0);

        // Check missing default
        assert_eq!(m.get(0, 2), -1.0);
        assert_eq!(m.get(3, 4), -1.0);
    }

    #[test]
    fn test_named_matrix_basic() {
        let names = vec!["A".to_string(), "B".to_string()];
        let mut m = NamedMatrix::new(names);
        
        m.set(0, 1, 0.5);
        
        // Check size
        assert_eq!(m.size(), 2);
        
        // Check values
        assert_eq!(m.get(0, 1), 0.5);
        assert_eq!(m.get(1, 0), 0.5); // Symmetric update
        assert_eq!(m.get(0, 0), 0.0);
        
        // Check by name
        assert_eq!(m.get_by_name("A", "B"), Some(0.5));
        assert_eq!(m.get_by_name("B", "A"), Some(0.5));
        assert_eq!(m.get_by_name("A", "C"), None);
        
        // Set by name
        assert!(m.set_by_name("A", "B", 0.8).is_ok());
        assert_eq!(m.get(0, 1), 0.8);
        assert!(m.set_by_name("A", "C", 0.9).is_err());
    }
}
