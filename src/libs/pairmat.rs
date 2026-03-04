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

/// A named matrix for storing pairwise distances/scores with sequence names.
///
/// Wraps a `CondensedMatrix` internally to save memory (N(N-1)/2).
/// Assumes symmetric matrix with 0 diagonal (distance matrix).
#[derive(Debug, Clone)]
pub struct NamedMatrix {
    names: indexmap::IndexMap<String, usize>,
    matrix: CondensedMatrix,
}

impl NamedMatrix {
    pub fn new(names: Vec<String>) -> Self {
        let size = names.len();
        let matrix = CondensedMatrix::new(size);
        let names: indexmap::IndexMap<_, _> = names
            .into_iter()
            .enumerate()
            .map(|(i, name)| (name, i))
            .collect();

        NamedMatrix { names, matrix }
    }

    /// Create from existing names and values (condensed upper triangle).
    pub fn new_from_values(names: Vec<String>, values: Vec<f32>) -> Self {
        let size = names.len();
        let matrix = CondensedMatrix::from_vec(size, values);
        
        let names: indexmap::IndexMap<_, _> = names
            .into_iter()
            .enumerate()
            .map(|(i, name)| (name, i))
            .collect();

        NamedMatrix { names, matrix }
    }

    /// Create with numeric names ("0", "1", ...).
    pub fn with_ids(size: usize) -> Self {
        let names: Vec<String> = (0..size).map(|i| i.to_string()).collect();
        Self::new(names)
    }

    pub fn size(&self) -> usize {
        self.matrix.size()
    }

    /// Access the underlying CondensedMatrix
    pub fn matrix(&self) -> &CondensedMatrix {
        &self.matrix
    }

    pub fn get(&self, row: usize, col: usize) -> f32 {
        self.matrix.get(row, col)
    }

    pub fn set(&mut self, row: usize, col: usize, value: f32) {
        self.matrix.set(row, col, value)
    }

    pub fn index(&self, row: usize, col: usize) -> usize {
        let (r, c) = if row < col { (row, col) } else { (col, row) };
        get_condensed_index(self.size(), r, c)
    }

    pub fn get_names(&self) -> Vec<&String> {
        self.names.keys().collect()
    }

    /// Get the underlying condensed data vector.
    pub fn values(&self) -> &[f32] {
        self.matrix.data()
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
        let mut matrix = NamedMatrix::new(index_name.into_iter().collect());
        for i in 0..size {
            for j in (i + 1)..size {
                matrix.set(i, j, scoring_matrix.get(i, j));
            }
        }
        matrix
    }

    /// Creates a new matrix from a relaxed PHYLIP format file
    ///
    /// ```no_run
    /// # use pgr::libs::pairmat::NamedMatrix;
    /// let matrix = NamedMatrix::from_relaxed_phylip("input.phy");
    /// ```
    pub fn from_relaxed_phylip(infile: &str) -> Self {
        let mut names = Vec::new();
        let mut raw_values = Vec::new();

        let reader = crate::reader(infile);
        let mut lines = reader.lines();

        // Skip the optional sequence count line
        if let Some(Ok(line)) = lines.next() {
            if line.trim().parse::<usize>().is_err() {
                // If first line is not a number, treat it as a data line
                Self::process_phylip_line(&line, &mut names, &mut raw_values);
            }
        }

        // Process remaining lines
        for line in lines.map_while(Result::ok) {
            Self::process_phylip_line(&line, &mut names, &mut raw_values);
        }

        let size = names.len();
        let mut matrix = Self::new(names);

        // Fill the matrix (lower triangle from PHYLIP)
        // raw_values contains flattened lower triangle: (1,0), (2,0), (2,1), ...
        let mut k = 0;
        for i in 0..size {
            for j in 0..=i {
                if k < raw_values.len() {
                    let value = raw_values[k];
                    matrix.set(i, j, value);
                    k += 1;
                }
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

/// Convert row, col (where row < col) to linear index in condensed array (upper triangle).
///
/// Based on formula: k = N*row - row*(row+1)/2 + col - row - 1
#[inline]
pub fn get_condensed_index(size: usize, row: usize, col: usize) -> usize {
    debug_assert!(row < col);
    debug_assert!(col < size);
    size * row - (row * (row + 1)) / 2 + col - row - 1
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
        assert_eq!(get_condensed_index(m.size(), 0, 1), 0);
        assert_eq!(get_condensed_index(m.size(), 0, 2), 1);
        assert_eq!(get_condensed_index(m.size(), 0, 3), 2);
        assert_eq!(get_condensed_index(m.size(), 1, 2), 3);
        assert_eq!(get_condensed_index(m.size(), 1, 3), 4);
        assert_eq!(get_condensed_index(m.size(), 2, 3), 5);
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
        
        // Test underlying data access
        let data = m.data();
        assert_eq!(data.len(), 3); // 3*2/2 = 3
        // Order: (0,1), (0,2), (1,2) -> 1.0, 2.0, 3.0
        assert_eq!(data[0], 1.0);
        assert_eq!(data[1], 2.0);
        assert_eq!(data[2], 3.0);
    }

    #[test]
    fn test_condensed_matrix_from_vec() {
        let data = vec![1.0, 2.0, 3.0];
        let m = CondensedMatrix::from_vec(3, data);
        assert_eq!(m.get(0, 1), 1.0);
        assert_eq!(m.get(0, 2), 2.0);
        assert_eq!(m.get(1, 2), 3.0);
    }

    #[test]
    #[should_panic(expected = "Data length 2 does not match expected length 3 for size 3")]
    fn test_condensed_matrix_from_vec_invalid_len() {
        CondensedMatrix::from_vec(3, vec![1.0, 2.0]);
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
        assert_eq!(m.get(0, 1), 0.5);
        assert_eq!(m.get(1, 0), 0.5);
        assert_eq!(m.get(0, 0), 0.0);
        
        assert_eq!(m.get_by_name("A", "B"), Some(0.5));
    }
    
    #[test]
    fn test_named_matrix_indexing() {
        let names = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        let m = NamedMatrix::new(names);
        
        // Size 3 -> len 3
        assert_eq!(m.values().len(), 3);
        
        // Index check
        // (0,1) -> 0
        // (0,2) -> 1
        // (1,2) -> 2
        assert_eq!(m.index(0, 1), 0);
        assert_eq!(m.index(0, 2), 1);
        assert_eq!(m.index(1, 2), 2);
    }
}
