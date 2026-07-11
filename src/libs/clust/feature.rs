//! Feature vector for distance-based clustering.
//!
//! A `FeatureVector` pairs a name with a list of float coordinates, used
//! by `cmd_pgr/dist/vector.rs` and `libs/clust/eval.rs` for distance
//! computation and cluster evaluation.

use anyhow::anyhow;
use std::io::BufRead;

//----------------------------
// FeatureVector
//----------------------------
#[derive(Default, Clone)]
pub struct FeatureVector {
    name: String,
    list: Vec<f32>,
}

impl FeatureVector {
    // Immutable accessors
    pub fn name(&self) -> &String {
        &self.name
    }
    pub fn list(&self) -> &Vec<f32> {
        &self.list
    }

    pub fn new() -> Self {
        Self {
            name: String::new(),
            list: vec![],
        }
    }

    /// Constructed from range and seq
    ///
    /// ```ignore
    /// # use pgr::libs::clust::feature::FeatureVector;
    /// let name = "Es_coli_005008_GCF_013426115_1".to_string();
    /// let list : Vec<f32> = vec![1.0,5.0,2.0,7.0,6.0,6.0];
    /// let entry = FeatureVector::from(&name, &list);
    /// # assert_eq!(*entry.name(), "Es_coli_005008_GCF_013426115_1");
    /// # assert_eq!(*entry.list().get(1).unwrap(), 5f32);
    /// ```
    pub fn from(name: &str, vector: &[f32]) -> Self {
        Self {
            name: name.to_owned(),
            list: Vec::from(vector),
        }
    }

    /// ```ignore
    /// # use pgr::libs::clust::feature::FeatureVector;
    /// let line = "Es_coli_005008_GCF_013426115_1\t1,5,2,7,6,6".to_string();
    /// let entry = FeatureVector::parse(&line).unwrap();
    /// # assert_eq!(*entry.name(), "Es_coli_005008_GCF_013426115_1");
    /// # assert_eq!(*entry.list().get(1).unwrap(), 5f32);
    /// ```
    pub fn parse(line: &str) -> anyhow::Result<FeatureVector> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() == 2 {
            let name = fields[0].to_string();
            let parts: Vec<&str> = fields[1].split(',').collect();
            let list: Vec<f32> = parts
                .iter()
                .map(|e| {
                    e.parse::<f32>()
                        .map_err(|e| anyhow!("invalid float value: {}", e))
                })
                .collect::<anyhow::Result<_>>()?;
            Ok(Self::from(&name, &list))
        } else {
            Ok(Self::new())
        }
    }
}

impl std::fmt::Display for FeatureVector {
    /// To string
    ///
    /// ```ignore
    /// # use pgr::libs::clust::feature::FeatureVector;
    /// let name = "Es_coli_005008_GCF_013426115_1".to_string();
    /// let list : Vec<f32> = vec![1.0,5.0,2.0,7.0,6.0,6.0];
    /// let entry = FeatureVector::from(&name, &list);
    /// assert_eq!(entry.to_string(), "Es_coli_005008_GCF_013426115_1\t1,5,2,7,6,6\n");
    /// ```
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(
            f,
            "{}\t{}",
            self.name(),
            self.list
                .iter()
                .map(|e| e.to_string())
                .collect::<Vec<_>>()
                .join(","),
        )?;
        Ok(())
    }
}

/// Load feature vectors from a file, optionally binarizing values to 0.0/1.0.
pub fn load_feature_vectors(infile: &str, is_bin: bool) -> anyhow::Result<Vec<FeatureVector>> {
    let mut entries = vec![];
    let reader = crate::libs::io::reader(infile)?;
    for line in reader.lines() {
        let line = line?;
        let mut entry = FeatureVector::parse(&line)?;
        if entry.name().is_empty() {
            continue;
        }
        if is_bin {
            let bin_list = entry
                .list()
                .iter()
                .map(|e| if *e > 0.0 { 1.0 } else { 0.0 })
                .collect::<Vec<f32>>();
            entry = FeatureVector::from(entry.name(), &bin_list);
        }
        entries.push(entry);
    }
    Ok(entries)
}
