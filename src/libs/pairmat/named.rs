use std::io::BufRead;

use super::condensed::{get_condensed_index, CondensedMatrix};
use super::scoring::ScoringMatrix;

/// A named matrix for storing pairwise distances/scores with sequence names.
///
/// Wraps a `CondensedMatrix` internally to save memory (N(N-1)/2).
/// Assumes symmetric matrix with 0 diagonal (distance matrix).
#[derive(Debug, Clone)]
pub struct NamedMatrix {
    names: indexmap::IndexMap<String, usize>,
    matrix: CondensedMatrix,
    diags: Option<Vec<f32>>,
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

        NamedMatrix {
            names,
            matrix,
            diags: None,
        }
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

        NamedMatrix {
            names,
            matrix,
            diags: None,
        }
    }

    /// Create with numeric names ("0", "1", ...).
    pub fn with_ids(size: usize) -> Self {
        let names: Vec<String> = (0..size).map(|i| i.to_string()).collect();
        Self::new(names)
    }

    pub fn size(&self) -> usize {
        self.matrix.size()
    }

    /// Consume the NamedMatrix and return its parts (names, condensed matrix).
    pub fn into_parts(self) -> (Vec<String>, CondensedMatrix) {
        let names = self.names.into_keys().collect();
        (names, self.matrix)
    }

    /// Access the underlying CondensedMatrix
    pub fn matrix(&self) -> &CondensedMatrix {
        &self.matrix
    }

    pub fn get(&self, row: usize, col: usize) -> f32 {
        if row == col {
            if let Some(ref diags) = self.diags {
                return diags[row];
            }
        }
        self.matrix.get(row, col)
    }

    pub fn set(&mut self, row: usize, col: usize, value: f32) {
        if row == col {
            if let Some(ref mut diags) = self.diags {
                diags[row] = value;
            }
        } else {
            self.matrix.set(row, col, value)
        }
    }

    pub fn index(&self, row: usize, col: usize) -> usize {
        let (r, c) = if row < col { (row, col) } else { (col, row) };
        get_condensed_index(self.size(), r, c)
    }

    pub fn get_names(&self) -> Vec<&String> {
        self.names.keys().collect()
    }

    pub fn get_index(&self, name: &str) -> Option<usize> {
        self.names.get(name).copied()
    }

    pub fn set_diags(&mut self, diags: Vec<f32>) {
        if diags.len() == self.size() {
            self.diags = Some(diags);
        }
    }

    pub fn get_diags(&self) -> Option<&Vec<f32>> {
        self.diags.as_ref()
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

    pub fn from_pair_scores(infile: &str, same: f32, missing: f32) -> anyhow::Result<Self> {
        let (scoring_matrix, index_name) = ScoringMatrix::from_pair_scores(infile, same, missing)?;
        let size = index_name.len();

        // Create NamedMatrix from ScoringMatrix
        let mut matrix = NamedMatrix::new(index_name.into_iter().collect());
        let mut diags = vec![same; size];

        for (i, d) in diags.iter_mut().enumerate().take(size) {
            *d = scoring_matrix.get(i, i);
            for j in (i + 1)..size {
                matrix.set(i, j, scoring_matrix.get(i, j));
            }
        }
        matrix.set_diags(diags);
        Ok(matrix)
    }

    /// Creates a new matrix from a relaxed PHYLIP format file
    ///
    /// ```no_run
    /// # use pgr::libs::pairmat::NamedMatrix;
    /// let matrix = NamedMatrix::from_relaxed_phylip("input.phy").unwrap();
    /// ```
    pub fn from_relaxed_phylip(infile: &str) -> anyhow::Result<Self> {
        let mut names = Vec::new();
        let mut raw_values = Vec::new();

        let reader = crate::reader(infile)?;
        let mut lines = reader.lines();

        // Skip the optional sequence count line
        if let Some(Ok(line)) = lines.next() {
            if line.trim().parse::<usize>().is_err() {
                // If first line is not a number, treat it as a data line
                Self::process_phylip_line(&line, &mut names, &mut raw_values)?;
            }
        }

        // Process remaining lines
        for line in lines.map_while(Result::ok) {
            Self::process_phylip_line(&line, &mut names, &mut raw_values)?;
        }

        let size = names.len();
        let mut matrix = Self::new(names);
        let mut diags = vec![0.0; size];

        // Fill the matrix (lower triangle from PHYLIP)
        // raw_values contains flattened lower triangle: (1,0), (2,0), (2,1), ...
        let mut k = 0;
        for (i, d) in diags.iter_mut().enumerate().take(size) {
            for j in 0..=i {
                if k < raw_values.len() {
                    let value = raw_values[k];
                    if i == j {
                        *d = value;
                    } else {
                        matrix.set(i, j, value);
                    }
                    k += 1;
                }
            }
        }
        matrix.set_diags(diags);
        Ok(matrix)
    }

    fn process_phylip_line(
        line: &str,
        names: &mut Vec<String>,
        values: &mut Vec<f32>,
    ) -> anyhow::Result<()> {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if !parts.is_empty() {
            let name = parts[0].to_string();
            names.push(name);

            // Read lower-triangle distances
            let distances: Vec<f32> = parts[1..=names.len()]
                .iter()
                .map(|&s| {
                    s.parse::<f32>()
                        .map_err(|e| anyhow::anyhow!("parse error: {e}"))
                })
                .collect::<anyhow::Result<Vec<f32>>>()?;

            values.extend(distances);
        }
        Ok(())
    }
}
