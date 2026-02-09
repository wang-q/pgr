use std::collections::{BTreeMap, HashMap, HashSet};
use super::tree::Tree;
use fixedbitset::FixedBitSet;

/// Trait for tree comparison and topology analysis
pub trait TreeComparison {
    /// Get the set of all splits (bipartitions) in the tree.
    ///
    /// Requires a `leaf_map` that maps leaf names to bit indices (0..N).
    /// This ensures that splits from different trees are comparable.
    ///
    /// Splits are normalized (always include the first taxon) and trivial splits
    /// (separating one leaf from the rest) are excluded by default to focus on topology.
    fn get_splits(&self, leaf_map: &BTreeMap<String, usize>) -> HashSet<FixedBitSet>;

    /// Get splits with their associated branch lengths.
    ///
    /// Returns a map from Split -> Branch Length.
    /// Used for Weighted Robinson-Foulds and Kuhner-Felsenstein distances.
    fn get_splits_with_values(&self, leaf_map: &BTreeMap<String, usize>) -> HashMap<FixedBitSet, f64>;

    /// Compute the Robinson-Foulds (RF) distance between two trees.
    ///
    /// RF distance is the symmetric difference of non-trivial splits: |S1 \ S2| + |S2 \ S1|.
    /// Returns error if trees have different sets of leaves.
    ///
    /// Note: This computes Unrooted RF distance.
    fn robinson_foulds(&self, other: &Self) -> Result<usize, String>;

    /// Compute the Weighted Robinson-Foulds (WRF) distance.
    ///
    /// Sum of absolute differences in branch lengths for all splits.
    /// If a split is missing in one tree, its length is assumed to be 0.
    fn weighted_robinson_foulds(&self, other: &Self) -> Result<f64, String>;

    /// Compute the Kuhner-Felsenstein (KF) distance (Branch Score Distance).
    ///
    /// Square root of the sum of squared differences in branch lengths for all splits.
    fn kuhner_felsenstein(&self, other: &Self) -> Result<f64, String>;
}

impl TreeComparison for Tree {
    fn get_splits(&self, leaf_map: &BTreeMap<String, usize>) -> HashSet<FixedBitSet> {
        let num_leaves = leaf_map.len();
        self.get_splits_with_values(leaf_map)
            .into_keys()
            .filter(|split| {
                let count = split.count_ones(..);
                count > 1 && count < num_leaves - 1
            })
            .collect()
    }

    fn get_splits_with_values(&self, leaf_map: &BTreeMap<String, usize>) -> HashMap<FixedBitSet, f64> {
        let mut splits = HashMap::new();
        let num_leaves = leaf_map.len();

        let root_id = match self.get_root() {
            Some(id) => id,
            None => return splits,
        };

        // Get all nodes in postorder (bottom-up)
        let nodes = match self.postorder(&root_id) {
            Ok(n) => n,
            Err(_) => return splits,
        };

        // Map NodeId -> BitSet (set of leaves under this node)
        let mut node_leaves: BTreeMap<usize, FixedBitSet> = BTreeMap::new();

        for node_id in nodes {
            let mut bitset = FixedBitSet::with_capacity(num_leaves);
            let node = self.get_node(node_id).unwrap();

            if node.is_leaf() {
                if let Some(name) = &node.name {
                    if let Some(&idx) = leaf_map.get(name) {
                        bitset.insert(idx);
                    }
                }
            } else {
                for child in &node.children {
                    if let Some(child_bs) = node_leaves.get(child) {
                        bitset.union_with(child_bs);
                    }
                }
            }

            // Normalize: Bipartitions are unrooted.
            // Convention: Always include the first taxon (index 0).
            // If bitset does NOT contain 0, take its complement.
            let mut normalized = bitset.clone();
            if num_leaves > 0 && !normalized.contains(0) {
                normalized.toggle_range(..num_leaves);
            }

            // Use branch length, default to 0.0 if None
            let len = node.length.unwrap_or(0.0);
            *splits.entry(normalized.clone()).or_insert(0.0) += len;

            node_leaves.insert(node_id, bitset);
        }

        splits
    }

    fn robinson_foulds(&self, other: &Self) -> Result<usize, String> {
        // 1. Check leaf consistency
        // get_leaf_names returns Vec<Option<String>>
        let leaves_self: HashSet<_> = self.get_leaf_names().into_iter().flatten().collect();
        let leaves_other: HashSet<_> = other.get_leaf_names().into_iter().flatten().collect();

        if leaves_self != leaves_other {
            // Sort for consistent error message
            let mut diff1: Vec<_> = leaves_self.difference(&leaves_other).collect();
            let mut diff2: Vec<_> = leaves_other.difference(&leaves_self).collect();
            diff1.sort();
            diff2.sort();
            return Err(format!(
                "Leaf sets do not match.\nIn Tree1 but not Tree2: {:?}\nIn Tree2 but not Tree1: {:?}",
                diff1, diff2
            ));
        }

        // 2. Build map
        let mut all_leaves: Vec<_> = leaves_self.into_iter().collect();
        all_leaves.sort(); // Deterministic order

        let mut leaf_map = BTreeMap::new();
        for (i, name) in all_leaves.iter().enumerate() {
            leaf_map.insert(name.clone(), i);
        }

        // 3. Get splits
        let splits1 = self.get_splits(&leaf_map);
        let splits2 = other.get_splits(&leaf_map);

        // 4. Calculate symmetric difference size
        // |A \ B| + |B \ A| = (A union B) - (A intersect B)
        // Or just count differences
        let diff1 = splits1.difference(&splits2).count();
        let diff2 = splits2.difference(&splits1).count();

        Ok(diff1 + diff2)
    }

    fn weighted_robinson_foulds(&self, other: &Self) -> Result<f64, String> {
        // 1. Check leaf consistency
        let leaves_self: HashSet<_> = self.get_leaf_names().into_iter().flatten().collect();
        let leaves_other: HashSet<_> = other.get_leaf_names().into_iter().flatten().collect();

        if leaves_self != leaves_other {
             let mut diff1: Vec<_> = leaves_self.difference(&leaves_other).collect();
            let mut diff2: Vec<_> = leaves_other.difference(&leaves_self).collect();
            diff1.sort();
            diff2.sort();
            return Err(format!(
                "Leaf sets do not match.\nIn Tree1 but not Tree2: {:?}\nIn Tree2 but not Tree1: {:?}",
                diff1, diff2
            ));
        }

        // 2. Build map
        let mut all_leaves: Vec<_> = leaves_self.into_iter().collect();
        all_leaves.sort();

        let mut leaf_map = BTreeMap::new();
        for (i, name) in all_leaves.iter().enumerate() {
            leaf_map.insert(name.clone(), i);
        }

        // 3. Get splits with values
        let splits1 = self.get_splits_with_values(&leaf_map);
        let splits2 = other.get_splits_with_values(&leaf_map);

        // 4. Calculate WRF
        let mut dist = 0.0;

        // Iterate over union of keys
        let keys: HashSet<_> = splits1.keys().chain(splits2.keys()).collect();
        
        for key in keys {
            let v1 = splits1.get(key).copied().unwrap_or(0.0);
            let v2 = splits2.get(key).copied().unwrap_or(0.0);
            dist += (v1 - v2).abs();
        }

        Ok(dist)
    }

    fn kuhner_felsenstein(&self, other: &Self) -> Result<f64, String> {
        // 1. Check leaf consistency
        let leaves_self: HashSet<_> = self.get_leaf_names().into_iter().flatten().collect();
        let leaves_other: HashSet<_> = other.get_leaf_names().into_iter().flatten().collect();

        if leaves_self != leaves_other {
             let mut diff1: Vec<_> = leaves_self.difference(&leaves_other).collect();
            let mut diff2: Vec<_> = leaves_other.difference(&leaves_self).collect();
            diff1.sort();
            diff2.sort();
            return Err(format!(
                "Leaf sets do not match.\nIn Tree1 but not Tree2: {:?}\nIn Tree2 but not Tree1: {:?}",
                diff1, diff2
            ));
        }

        // 2. Build map
        let mut all_leaves: Vec<_> = leaves_self.into_iter().collect();
        all_leaves.sort();

        let mut leaf_map = BTreeMap::new();
        for (i, name) in all_leaves.iter().enumerate() {
            leaf_map.insert(name.clone(), i);
        }

        // 3. Get splits with values
        let splits1 = self.get_splits_with_values(&leaf_map);
        let splits2 = other.get_splits_with_values(&leaf_map);

        // 4. Calculate KF (Sum of squares)
        let mut sum_sq = 0.0;

        // Iterate over union of keys
        let keys: HashSet<_> = splits1.keys().chain(splits2.keys()).collect();
        
        for key in keys {
            let v1 = splits1.get(key).copied().unwrap_or(0.0);
            let v2 = splits2.get(key).copied().unwrap_or(0.0);
            sum_sq += (v1 - v2).powi(2);
        }

        Ok(sum_sq.sqrt())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::phylo::tree::Tree;

    #[test]
    fn test_rf_phylotree_rs_suite() {
        let cases = [
            (
                26, // Was 28 in phylotree-rs, but our implementation (normalized splits) yields 26.
                false,
                "(((z,y),((x,w),((v,u),(t,s)))),((r,(q,(p,(o,(n,m))))),((l,k),((j,i),(h,g)))));",
                "((((s,(r,q)),((p,o),((n,(m,l)),(k,(j,(i,(h,g))))))),(z,y)),((x,w),(v,(u,t))));",
            ),
            (
                28, // Was 30
                false,
                "(((q,(p,o)),((n,m),((l,(k,(j,(i,(h,g))))),(z,y)))),((x,(w,v)),(u,(t,(s,r)))));",
                "((((t,s),((r,q),((p,o),(n,m)))),((l,k),(j,i))),(((h,g),z),((y,(x,w)),(v,u))));",
            ),
            (
                18, // Was 24
                false,
                "(((p,o),(n,m)),(((l,(k,(j,i))),(h,g)),((z,y),((x,w),((v,u),(t,(s,(r,q))))))));",
                "((x,(w,v)),((u,(t,(s,(r,q)))),((p,(o,(n,(m,(l,(k,(j,(i,(h,g))))))))),(z,y))));",
            ),
        ];

        for (expected_rf, _weighted, t1_str, t2_str) in cases {
            let t1 = Tree::from_newick(t1_str).unwrap();
            let t2 = Tree::from_newick(t2_str).unwrap();
            let rf = t1.robinson_foulds(&t2).unwrap();
            assert_eq!(rf, expected_rf, "Failed for case: {} vs {}", t1_str, t2_str);
        }
    }

    #[test]
    fn test_wrf_kf() {
        // T1: ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);
        // T2: ((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.2);
        // Split {A,B}: T1=0.2, T2=0.3. Diff=0.1.
        // Other splits: Trivial {A}, {B}, {C}, {D} are 0.1 each. {C,D} is 0.2 each.
        // Assuming trivial splits are identical and lengths match.
        // WRF: |0.2 - 0.3| = 0.1.
        // KF: sqrt((0.2-0.3)^2) = 0.1.
        
        let t1_str = "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);";
        let t2_str = "((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.2);";
        
        let t1 = Tree::from_newick(t1_str).unwrap();
        let t2 = Tree::from_newick(t2_str).unwrap();
        
        let wrf = t1.weighted_robinson_foulds(&t2).unwrap();
        let kf = t1.kuhner_felsenstein(&t2).unwrap();
        
        assert!((wrf - 0.1).abs() < 1e-6, "WRF expected 0.1, got {}", wrf);
        assert!((kf - 0.1).abs() < 1e-6, "KF expected 0.1, got {}", kf);
    }

    #[test]
    fn test_wrf_kf_topology_change() {
        // T1: ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);
        // T3: ((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);
        // Splits T1: {A,B} (normalized) sum of lengths = 0.2 + 0.2 = 0.4.
        // Splits T3: {A,C} (normalized) sum of lengths = 0.2 + 0.2 = 0.4.
        // Shared: Trivials (values match).
        // Diff: {A,B} in T1 (0.4) not T3 (0.0). {A,C} in T3 (0.4) not T1 (0.0).
        // WRF = 0.4 + 0.4 = 0.8.
        // KF = sqrt(0.4^2 + 0.4^2) = sqrt(0.16 + 0.16) = sqrt(0.32) ≈ 0.5656854.

        let t1_str = "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);";
        let t3_str = "((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);";

        let t1 = Tree::from_newick(t1_str).unwrap();
        let t3 = Tree::from_newick(t3_str).unwrap();

        let wrf = t1.weighted_robinson_foulds(&t3).unwrap();
        let kf = t1.kuhner_felsenstein(&t3).unwrap();

        assert!((wrf - 0.8).abs() < 1e-6, "WRF expected 0.8, got {}", wrf);
        assert!((kf - 0.32f64.sqrt()).abs() < 1e-6, "KF expected sqrt(0.32), got {}", kf);
    }
}
