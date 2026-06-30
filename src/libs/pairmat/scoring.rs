use std::collections::HashMap;
use std::io::BufRead;

/// A symmetric scoring matrix parameterized over the value type `T`.
#[derive(Debug, Clone)]
pub struct ScoringMatrix<T> {
    size: Option<usize>,
    same: Option<T>,
    missing: Option<T>,
    data: HashMap<(usize, usize), T>,
}

impl<T> Default for ScoringMatrix<T>
where
    T: Default + Copy,
{
    fn default() -> Self {
        Self::new()
    }
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
                .unwrap_or_else(|| self.same.unwrap_or_default())
        } else {
            let (r, c) = if row < col { (row, col) } else { (col, row) };
            self.data
                .get(&(r, c))
                .copied()
                .unwrap_or_else(|| self.missing.unwrap_or_default())
        }
    }
}

// Add a separate implementation for f32 specifically for from_pair_scores
impl ScoringMatrix<f32> {
    pub fn from_pair_scores(
        infile: &str,
        same: f32,
        missing: f32,
    ) -> anyhow::Result<(Self, Vec<String>)> {
        let mut names = indexmap::IndexSet::new();
        let mut matrix = Self::with_defaults(same, missing);

        let reader = crate::reader(infile)?;
        for line in reader.lines().map_while(Result::ok) {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 3 {
                let n1 = fields[0].to_string();
                let n2 = fields[1].to_string();
                let score: f32 = fields[2].parse()?;

                names.insert(n1.clone());
                names.insert(n2.clone());

                let i1 = names
                    .get_index_of(&n1)
                    .ok_or_else(|| anyhow::anyhow!("name not found: {n1}"))?;
                let i2 = names
                    .get_index_of(&n2)
                    .ok_or_else(|| anyhow::anyhow!("name not found: {n2}"))?;
                matrix.set(i1, i2, score);
            }
        }

        matrix.set_size(names.len());
        Ok((matrix, names.into_iter().collect()))
    }
}
