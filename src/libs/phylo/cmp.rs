
use super::tree::Tree;
use fixedbitset::FixedBitSet;
use std::collections::{BTreeMap, HashSet};

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

    /// Compute the Robinson-Foulds (RF) distance between two trees.
    ///
    /// RF distance is the symmetric difference of non-trivial splits: |S1 \ S2| + |S2 \ S1|.
    /// Returns error if trees have different sets of leaves.
    ///
    /// Note: This computes Unrooted RF distance.
    fn robinson_foulds(&self, other: &Self) -> Result<usize, String>;
}

impl TreeComparison for Tree {
    fn get_splits(&self, leaf_map: &BTreeMap<String, usize>) -> HashSet<FixedBitSet> {
        let mut splits = HashSet::new();
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

            // Filter trivial splits (optional, but standard for RF distance)
            // Trivial = size 0, size 1, size N, size N-1 (since we normalize, checks are simpler)
            // After normalization, 0 is always present.
            // So we only have sets containing 0.
            // Trivial cases containing 0:
            // - Size 1 (only 0) -> Leaf 0.
            // - Size N (all) -> Root.
            // - Size N-1 (all except one) -> Complement of some other leaf.
            
            let count = normalized.count_ones(..);
            let is_trivial = count <= 1 || count >= num_leaves - 1;

            if !is_trivial {
                splits.insert(normalized.clone());
            }

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
            diff1.sort();
            let mut diff2: Vec<_> = leaves_other.difference(&leaves_self).collect();
            diff2.sort();
            
            return Err(format!(
                "Trees have different leaf sets.\nIn Tree1 only: {:?}\nIn Tree2 only: {:?}",
                diff1, diff2
            ));
        }

        if leaves_self.is_empty() {
            return Ok(0);
        }

        // 2. Build canonical map
        // Sort leaves to ensure consistent indexing
        let mut sorted_leaves: Vec<_> = leaves_self.into_iter().collect();
        sorted_leaves.sort();

        let mut leaf_map = BTreeMap::new();
        for (i, name) in sorted_leaves.into_iter().enumerate() {
            leaf_map.insert(name, i);
        }

        // 3. Compute splits
        let splits_self = self.get_splits(&leaf_map);
        let splits_other = other.get_splits(&leaf_map);

        // 4. Symmetric Difference Size
        // RF = |S1 \ S2| + |S2 \ S1|
        let intersection_count = splits_self.intersection(&splits_other).count();
        let rf = splits_self.len() + splits_other.len() - 2 * intersection_count;

        Ok(rf)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::phylo::tree::Tree;

    #[test]
    fn test_rf_distance_identical() {
        let t1 = Tree::from_newick("((A,B),C);").unwrap();
        let t2 = Tree::from_newick("((A,B),C);").unwrap();
        assert_eq!(t1.robinson_foulds(&t2).unwrap(), 0);
    }

    #[test]
    fn test_rf_distance_different_topology() {
        let t1 = Tree::from_newick("((A,B),C);").unwrap(); 
        let t2 = Tree::from_newick("((A,C),B);").unwrap(); 
        // 3 leaves. Unrooted. No internal edges.
        // Expect 0.
        // Wait. ((A,B),C) has no internal edges.
        // ((A,C),B) has no internal edges.
        // So they are topologically identical as unrooted trees.
        assert_eq!(t1.robinson_foulds(&t2).unwrap(), 0);
    }

    #[test]
    fn test_rf_distance_star() {
        let t1 = Tree::from_newick("((A,B),C);").unwrap(); 
        let t2 = Tree::from_newick("(A,B,C);").unwrap();     
        // Both unrooted star trees.
        assert_eq!(t1.robinson_foulds(&t2).unwrap(), 0);
    }

    #[test]
    fn test_rf_distance_complex() {
        // Tree 1: ((A,B),(C,D));  
        // Internal Split: {A,B} vs {C,D}. Normalized to {A,B}.
        // Tree 2: ((A,C),(B,D));  
        // Internal Split: {A,C} vs {B,D}. Normalized to {A,C}.
        // Diff: {A,B} vs {A,C}.
        // RF = 2.
        let t1 = Tree::from_newick("((A,B),(C,D));").unwrap();
        let t2 = Tree::from_newick("((A,C),(B,D));").unwrap();
        assert_eq!(t1.robinson_foulds(&t2).unwrap(), 2);
    }
    
    #[test]
    fn test_rf_distance_5_taxa() {
        // (A,B,(C,D,E)) vs (A,B,(C,(D,E)))
        // T1: Star-like for C,D,E. No internal splits among C,D,E.
        // Internal splits: {A,B} vs {C,D,E}.
        // T2: {D,E} is a split. {C,D,E} is a split. {A,B} is a split.
        // T1 splits: {{A,B}} (normalized)
        // T2 splits: {{A,B}, {D,E}} (normalized)
        // Diff: {D,E} is in T2 not T1.
        // RF = 1.
        let t1 = Tree::from_newick("((A,B),(C,D,E));").unwrap();
        let t2 = Tree::from_newick("((A,B),(C,(D,E)));").unwrap();
        assert_eq!(t1.robinson_foulds(&t2).unwrap(), 1);
    }

    #[test]
    fn test_leaf_mismatch() {
        let t1 = Tree::from_newick("((A,B),C);").unwrap();
        let t2 = Tree::from_newick("((A,B),D);").unwrap();
        assert!(t1.robinson_foulds(&t2).is_err());
    }
}
