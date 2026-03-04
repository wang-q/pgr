use crate::libs::pairmat::{CondensedMatrix, NamedMatrix};
use crate::libs::phylo::tree::Tree;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Method {
    Single,
    Complete,
    Average,
    Weighted,
    Centroid,
    Median,
    Ward,
}

impl Method {
    /// Returns true if the method works in squared Euclidean space.
    /// These methods (Ward, Centroid, Median) require squared distances for linear updates.
    pub fn is_euclidean_squared(&self) -> bool {
        matches!(self, Method::Ward | Method::Centroid | Method::Median)
    }
}

impl std::str::FromStr for Method {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "single" => Ok(Method::Single),
            "complete" => Ok(Method::Complete),
            "average" | "upgma" => Ok(Method::Average),
            "weighted" | "wpgma" => Ok(Method::Weighted),
            "centroid" | "upgmc" => Ok(Method::Centroid),
            "median" | "wpgmc" => Ok(Method::Median),
            "ward" | "ward.d2" => Ok(Method::Ward),
            _ => Err(format!("Unknown linkage method: {}", s)),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Step {
    pub cluster1: usize,
    pub cluster2: usize,
    pub distance: f32,
    pub size: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Algorithm {
    Primitive,
    NnChain,
    Auto,
}

/// Perform hierarchical clustering on a condensed distance matrix.
///
/// Returns a list of steps (merges) forming the dendrogram.
/// The length of the result will be N-1 for a matrix of size N.
pub fn linkage(matrix: &NamedMatrix, method: Method) -> Vec<Step> {
    linkage_with_algo(matrix, method, Algorithm::Auto)
}

/// Perform hierarchical clustering with explicit algorithm selection.
pub fn linkage_with_algo(matrix: &NamedMatrix, method: Method, algo: Algorithm) -> Vec<Step> {
    // Create a mutable copy of the condensed matrix
    let mut condensed = CondensedMatrix::from_vec(matrix.size(), matrix.values().to_vec());

    linkage_core(&mut condensed, method, algo)
}

/// Perform hierarchical clustering in-place (consuming the distance matrix).
///
/// This avoids cloning the distance matrix, saving memory.
/// The input matrix will be modified (distances updated during merging).
pub fn linkage_inplace(mut matrix: CondensedMatrix, method: Method) -> Vec<Step> {
    linkage_core(&mut matrix, method, Algorithm::Auto)
}

fn linkage_core(condensed: &mut CondensedMatrix, method: Method, algo: Algorithm) -> Vec<Step> {
    match algo {
        Algorithm::Primitive => linkage_primitive(condensed, method),
        Algorithm::NnChain => linkage_nn_chain(condensed, method),
        Algorithm::Auto => match method {
            // Primitive O(N^3) is safest for Centroid/Median which don't satisfy reducibility
            Method::Centroid | Method::Median => linkage_primitive(condensed, method),

            // NN-chain O(N^2) for reducible metrics
            _ => linkage_nn_chain(condensed, method),
        },
    }
}

/// Nearest-neighbor chain algorithm for hierarchical clustering.
///
/// Time complexity: O(N^2)
/// Space complexity: O(N^2) (for distance matrix copy) + O(N) (chain stack)
///
/// Only valid for methods satisfying the reducibility property:
/// Single, Complete, Average, Weighted, Ward.
/// (Centroid and Median are NOT reducible and may cause infinite loops or incorrect results with NN-chain).
fn linkage_nn_chain(condensed: &mut CondensedMatrix, method: Method) -> Vec<Step> {
    let n = condensed.size();
    if n < 2 {
        return vec![];
    }

    // Optimization: Square the distances for Ward/Centroid/Median
    // This allows O(1) Lance-Williams updates without sqrt calls.
    let is_squared = method.is_euclidean_squared();
    if is_squared {
        for x in condensed.data_mut() {
            *x = *x * *x;
        }
    }

    // Cluster sizes
    let mut size = vec![1; n];

    // Active clusters (true if valid, false if merged)
    let mut active = vec![true; n];

    // Map internal index to original cluster ID (for step output)
    let mut cluster_ids: Vec<usize> = (0..n).collect();

    // The nearest-neighbor chain (stack of cluster indices)
    let mut chain = Vec::with_capacity(n);

    // Result steps
    let mut steps = Vec::with_capacity(n - 1);

    // Main loop
    // We need to perform N-1 merges.
    // Each merge reduces the number of active clusters by 1.
    // Instead of counting loop, we run until steps.len() == n - 1.
    while steps.len() < n - 1 {
        // If chain is empty, pick an arbitrary active cluster to start
        if chain.is_empty() {
            for i in 0..n {
                if active[i] {
                    chain.push(i);
                    break;
                }
            }
        }

        // Extend the chain from the top element
        // Find NN of chain.last()
        let k = *chain.last().unwrap();
        let mut min_dist = f32::INFINITY;
        let mut nn = k; // Default to self if no other active found (shouldn't happen if >1 active)

        for i in 0..n {
            if i == k || !active[i] {
                continue;
            }
            let d = condensed.get(k, i);
            // Strict inequality < is crucial for chain stability,
            // but for equal distances we need a tie-breaking rule (usually index).
            // Here: d < min_dist OR (d == min_dist and i < nn)
            if d < min_dist {
                min_dist = d;
                nn = i;
            }
        }

        // Check if NN(k) is already the previous element in chain (Reciprocal NN)
        if chain.len() >= 2 && nn == chain[chain.len() - 2] {
            // RNN found: merge k and nn
            let u = chain.pop().unwrap(); // k
            let v = chain.pop().unwrap(); // nn (which is also NN(k))

            // Ensure u < v for consistent indexing updates if needed,
            // though condensed matrix handles (u,v) order.
            // Let's stick to: we merge v into u (keep u, disable v) or vice versa.
            // Standard convention: keep the one with smaller index to minimize shifts?
            // Actually, we usually merge into one and mark other inactive.
            // Let's merge `u` and `v`.
            // Swap so u < v to keep smaller index active (arbitrary choice, but clean)
            let (u, v) = if u < v { (u, v) } else { (v, u) };

            // Distance between u and v
            let d_uv = min_dist;
            let dist_out = if is_squared { d_uv.sqrt() } else { d_uv };

            // Record step
            let id1 = cluster_ids[u];
            let id2 = cluster_ids[v];
            let size1 = size[u];
            let size2 = size[v];
            let new_size = size1 + size2;
            let new_id = n + steps.len(); // Next cluster ID

            steps.push(Step {
                cluster1: id1,
                cluster2: id2,
                distance: dist_out,
                size: new_size,
            });

            // Update distances (Lance-Williams)
            // Merge v into u
            for k in 0..n {
                if !active[k] || k == u || k == v {
                    continue;
                }

                let d_uk = condensed.get(u, k);
                let d_vk = condensed.get(v, k);

                let new_dist = lance_williams(method, d_uk, d_vk, d_uv, size1, size2, size[k]);

                condensed.set(u, k, new_dist);
            }

            // Update state
            size[u] = new_size;
            cluster_ids[u] = new_id;
            active[v] = false;

            // Since `v` is inactive, it must be removed from chain if present (it was just popped).
            // `u` is active and updated, it was also popped.
            // If `u` is still in chain (it isn't, we popped it), we'd need to check.
            // But we just popped both.
            // If `chain` is not empty, its new top might need to be re-evaluated against `u`.
            // So we loop back.
        } else {
            // Not RNN, push NN to chain
            chain.push(nn);
        }
    }

    // Post-processing: Sort steps by distance and re-assign IDs to match standard behavior
    // 1. Create indices and sort
    let mut indices: Vec<usize> = (0..steps.len()).collect();
    indices.sort_by(|&i, &j| {
        let s1 = &steps[i];
        let s2 = &steps[j];
        // Sort by distance ascending
        // If distances are equal, maintain original topological order (i vs j)
        s1.distance
            .partial_cmp(&s2.distance)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(i.cmp(&j))
    });

    // 2. Map old cluster IDs to new IDs
    // Leaves 0..N-1 are unchanged.
    // Internal nodes (steps) need remapping.
    // old_id = n + original_index
    // new_id = n + new_sorted_index
    let mut id_map: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for i in 0..n {
        id_map.insert(i, i);
    }

    let mut new_steps = Vec::with_capacity(steps.len());

    for (new_idx, &old_idx) in indices.iter().enumerate() {
        let step = &steps[old_idx];
        let old_res_id = n + old_idx;
        let new_res_id = n + new_idx;

        id_map.insert(old_res_id, new_res_id);

        let new_c1 = *id_map
            .get(&step.cluster1)
            .expect("Cluster ID not found in map");
        let new_c2 = *id_map
            .get(&step.cluster2)
            .expect("Cluster ID not found in map");

        // Ensure c1 < c2 for canonical output
        let (c1, c2) = if new_c1 < new_c2 {
            (new_c1, new_c2)
        } else {
            (new_c2, new_c1)
        };

        new_steps.push(Step {
            cluster1: c1,
            cluster2: c2,
            distance: step.distance,
            size: step.size,
        });
    }

    new_steps
}

/// Convert linkage steps to a Node (tree structure).
///
/// The resulting tree is rooted at the last merge.
/// Branch lengths are derived from merge heights.
pub fn to_tree(steps: &[Step], names: &[String]) -> Tree {
    let mut tree = Tree::new();
    let n = names.len();

    if n == 0 {
        return tree;
    }

    // Map original cluster ID to Tree NodeId
    // IDs 0..n-1 are leaves (original sequences)
    // IDs n..n+steps.len()-1 are internal nodes (merges)
    let mut cluster_to_node: std::collections::HashMap<usize, usize> =
        std::collections::HashMap::new();

    // Create leaf nodes
    for (i, name) in names.iter().enumerate() {
        let node_id = tree.add_node();
        if let Some(node) = tree.get_node_mut(node_id) {
            node.set_name(name.clone());
        }
        cluster_to_node.insert(i, node_id);
    }

    // Track height of each cluster (original ID -> height)
    let mut heights: std::collections::HashMap<usize, f32> = std::collections::HashMap::new();
    for i in 0..n {
        heights.insert(i, 0.0);
    }

    let mut root_cluster_id = 0;

    for (step_idx, step) in steps.iter().enumerate() {
        let new_cluster_id = n + step_idx;

        // Create parent node
        let parent_node_id = tree.add_node();
        cluster_to_node.insert(new_cluster_id, parent_node_id);

        // Process children
        for &child_cluster_id in &[step.cluster1, step.cluster2] {
            if let Some(&child_node_id) = cluster_to_node.get(&child_cluster_id) {
                let h_child = *heights.get(&child_cluster_id).unwrap_or(&0.0);

                // Calculate branch length
                // Use distance/2.0 as node height for ultrametric-like appearance
                let node_height = step.distance / 2.0;
                let len = (node_height - h_child).max(0.0); // Prevent negative length

                // Set length on child and link to parent
                if let Some(child_node) = tree.get_node_mut(child_node_id) {
                    child_node.length = Some(len as f64);
                }

                // Link in tree structure
                let _ = tree.add_child(parent_node_id, child_node_id);
            }
        }

        heights.insert(new_cluster_id, step.distance / 2.0);
        root_cluster_id = new_cluster_id;
    }

    // Set root of the tree
    if let Some(&root_node_id) = cluster_to_node.get(&root_cluster_id) {
        tree.set_root(root_node_id);
    } else if n == 1 {
        // Single leaf case
        if let Some(&root_node_id) = cluster_to_node.get(&0) {
            tree.set_root(root_node_id);
        }
    }

    tree
}

/// Primitive O(N^3) implementation of agglomerative clustering.
///
/// This serves as the MVP implementation (Phase 1) and a baseline for testing.
/// It maintains a mutable copy of the distance matrix and iteratively finds the minimum.
fn linkage_primitive(condensed: &mut CondensedMatrix, method: Method) -> Vec<Step> {
    let n = condensed.size();
    if n < 2 {
        return vec![];
    }

    // Optimization: Square the distances for Ward/Centroid/Median
    let is_squared = method.is_euclidean_squared();
    if is_squared {
        for x in condensed.data_mut() {
            *x = *x * *x;
        }
    }

    // Track cluster sizes and status
    let mut size = vec![1; n];
    let mut active = vec![true; n];

    // Map internal index to original cluster ID
    let mut cluster_ids: Vec<usize> = (0..n).collect();

    let mut steps = Vec::with_capacity(n - 1);

    for step_idx in 0..(n - 1) {
        // 1. Find min distance between active clusters
        let mut min_dist = f32::INFINITY;
        let mut u = 0;
        let mut v = 0;

        for i in 0..n {
            if !active[i] {
                continue;
            }
            for j in (i + 1)..n {
                if !active[j] {
                    continue;
                }
                let d = condensed.get(i, j);
                if d < min_dist {
                    min_dist = d;
                    u = i;
                    v = j;
                }
            }
        }

        // 2. Record the merge
        // New cluster ID
        let new_id = n + step_idx;
        let id1 = cluster_ids[u];
        let id2 = cluster_ids[v];
        let size1 = size[u];
        let size2 = size[v];
        let new_size = size1 + size2;

        let dist_out = if is_squared {
            min_dist.sqrt()
        } else {
            min_dist
        };

        steps.push(Step {
            cluster1: id1,
            cluster2: id2,
            distance: dist_out,
            size: new_size,
        });

        // 3. Update distances (Lance-Williams)
        // We merge v into u. u becomes the new cluster. v becomes inactive.
        for k in 0..n {
            if !active[k] || k == u || k == v {
                continue;
            }

            let d_uk = condensed.get(u, k);
            let d_vk = condensed.get(v, k);
            let d_uv = min_dist;

            let new_dist = lance_williams(method, d_uk, d_vk, d_uv, size1, size2, size[k]);

            condensed.set(u, k, new_dist);
        }

        // 4. Update state
        size[u] = new_size;
        cluster_ids[u] = new_id;
        active[v] = false;
    }

    steps
}

#[inline]
fn lance_williams(
    method: Method,
    d_uk: f32,
    d_vk: f32,
    d_uv: f32,
    size_u: usize,
    size_v: usize,
    size_k: usize,
) -> f32 {
    let n_u = size_u as f32;
    let n_v = size_v as f32;
    let n_k = size_k as f32;
    let n_uv = n_u + n_v;
    let n_uvk = n_uv + n_k;

    match method {
        Method::Single => d_uk.min(d_vk),
        Method::Complete => d_uk.max(d_vk),
        Method::Average => (n_u * d_uk + n_v * d_vk) / n_uv,
        Method::Weighted => 0.5 * (d_uk + d_vk),
        Method::Centroid => {
            // UPGMC (Unweighted Pair Group Method with Centroid Averaging)
            // Inputs are already squared distances.
            // d(u+v, k)^2 = (n_u * d(u,k)^2 + n_v * d(v,k)^2 - n_u*n_v*d(u,v)^2 / n_uv) / n_uv

            let d_new = (n_u * d_uk + n_v * d_vk - (n_u * n_v * d_uv) / n_uv) / n_uv;
            d_new.max(0.0)
        }
        Method::Median => {
            // WPGMC (Weighted Pair Group Method with Centroid Averaging)
            // Inputs are already squared distances.
            // d_new^2 = 0.5*d_uk^2 + 0.5*d_vk^2 - 0.25*d_uv^2
            let d_new = 0.5 * d_uk + 0.5 * d_vk - 0.25 * d_uv;
            d_new.max(0.0)
        }
        Method::Ward => {
            // Ward's method (minimal variance) - specifically Ward.D2
            // Inputs are already squared distances.
            // d_new^2 = ((n_u+n_k)*d_uk^2 + (n_v+n_k)*d_vk^2 - n_k*d_uv^2) / n_uvk

            let d_new = ((n_u + n_k) * d_uk + (n_v + n_k) * d_vk - n_k * d_uv) / n_uvk;
            d_new.max(0.0)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper to create a NamedMatrix for testing
    fn create_test_matrix(size: usize) -> NamedMatrix {
        let names: Vec<String> = (0..size).map(|i| i.to_string()).collect();
        NamedMatrix::new(names)
    }

    #[test]
    fn test_linkage_single() {
        // Distances:
        // 0-1: 1.0
        // 0-2: 4.0
        // 1-2: 2.0
        //
        // Steps:
        // 1. Merge 0-1 (d=1.0). New cluster 3 = {0,1}.
        //    Update dist(3, 2) = min(d(0,2), d(1,2)) = min(4, 2) = 2.0
        // 2. Merge 3-2 (d=2.0).

        let mut m = create_test_matrix(3);
        m.set(0, 1, 1.0);
        m.set(0, 2, 4.0);
        m.set(1, 2, 2.0);

        let steps = linkage(&m, Method::Single);

        assert_eq!(steps.len(), 2);

        assert_eq!(steps[0].cluster1, 0);
        assert_eq!(steps[0].cluster2, 1);
        assert_eq!(steps[0].distance, 1.0);
        assert_eq!(steps[0].size, 2);

        // Step 2: Merge {0,1} (id 3) and 2.
        // Canonical order: min(2, 3) = 2, max(2, 3) = 3.
        assert_eq!(steps[1].cluster1, 2);
        assert_eq!(steps[1].cluster2, 3); // 3 is the new id for {0,1}
        assert_eq!(steps[1].distance, 2.0);
        assert_eq!(steps[1].size, 3);
    }

    #[test]
    fn test_linkage_complete() {
        // Distances:
        // 0-1: 1.0
        // 0-2: 4.0
        // 1-2: 2.0
        //
        // Steps:
        // 1. Merge 0-1 (d=1.0). New cluster 3 = {0,1}.
        //    Update dist(3, 2) = max(d(0,2), d(1,2)) = max(4, 2) = 4.0
        // 2. Merge 3-2 (d=4.0).

        let mut m = create_test_matrix(3);
        m.set(0, 1, 1.0);
        m.set(0, 2, 4.0);
        m.set(1, 2, 2.0);

        let steps = linkage(&m, Method::Complete);

        assert_eq!(steps.len(), 2);
        assert_eq!(steps[1].distance, 4.0);
    }

    #[test]
    fn test_linkage_average() {
        // Distances:
        // 0-1: 1.0
        // 0-2: 4.0
        // 1-2: 2.0
        //
        // Steps:
        // 1. Merge 0-1 (d=1.0). Size=2.
        //    Update dist(3, 2) = (1*4.0 + 1*2.0) / 2 = 3.0
        // 2. Merge 3-2 (d=3.0).

        let mut m = create_test_matrix(3);
        m.set(0, 1, 1.0);
        m.set(0, 2, 4.0);
        m.set(1, 2, 2.0);

        let steps = linkage(&m, Method::Average);

        assert_eq!(steps.len(), 2);
        assert_eq!(steps[1].distance, 3.0);
    }

    #[test]
    fn test_nn_chain_vs_primitive() {
        // Create a random-ish matrix (5x5)
        // 0-1: 10
        // 0-2: 2
        // 0-3: 8
        // 0-4: 5
        // 1-2: 9
        // 1-3: 3
        // 1-4: 7
        // 2-3: 6
        // 2-4: 4
        // 3-4: 1

        let mut m = create_test_matrix(5);
        m.set(0, 1, 10.0);
        m.set(0, 2, 2.0);
        m.set(0, 3, 8.0);
        m.set(0, 4, 5.0);
        m.set(1, 2, 9.0);
        m.set(1, 3, 3.0);
        m.set(1, 4, 7.0);
        m.set(2, 3, 6.0);
        m.set(2, 4, 4.0);
        m.set(3, 4, 1.0); // min

        // Test with Average linkage (Reducible)
        let steps_prim = linkage_with_algo(&m, Method::Average, Algorithm::Primitive);
        let steps_nn = linkage_with_algo(&m, Method::Average, Algorithm::NnChain);

        assert_eq!(steps_prim.len(), 4);
        assert_eq!(steps_nn.len(), 4);

        for (i, (s1, s2)) in steps_prim.iter().zip(steps_nn.iter()).enumerate() {
            // Check distance (should be identical)
            assert!(
                (s1.distance - s2.distance).abs() < 1e-5,
                "Step {}: distance mismatch {} vs {}",
                i,
                s1.distance,
                s2.distance
            );

            // Check clusters (order might differ in representation but set should be same)
            // But for simple cases with strict inequality, they should be identical.
            // Let's check normalized cluster pairs
            let (min1, max1) = if s1.cluster1 < s1.cluster2 {
                (s1.cluster1, s1.cluster2)
            } else {
                (s1.cluster2, s1.cluster1)
            };
            let (min2, max2) = if s2.cluster1 < s2.cluster2 {
                (s2.cluster1, s2.cluster2)
            } else {
                (s2.cluster2, s2.cluster1)
            };

            assert_eq!(min1, min2, "Step {}: cluster1 mismatch", i);
            assert_eq!(max1, max2, "Step {}: cluster2 mismatch", i);
        }
    }

    // Helper to create a random matrix for testing
    fn create_random_matrix_local(size: usize) -> NamedMatrix {
        // Use a simple pseudo-random generator to avoid external deps or complexity
        let mut seed: u64 = 12345;
        let mut rng = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            (seed >> 33) as f32 / 2147483648.0
        };

        // Let's just create names and fill CondensedMatrix manually
        let names: Vec<String> = (0..size).map(|i| i.to_string()).collect();
        let mut m = NamedMatrix::new(names);

        for i in 0..size {
            for j in (i + 1)..size {
                m.set(i, j, rng());
            }
        }
        m
    }

    #[test]
    fn test_nn_chain_fuzzing() {
        // Run multiple random tests to ensure stability and correctness
        for i in 0..20 {
            let size = 10 + i * 5; // Sizes: 10, 15, ..., 105
            let m = create_random_matrix_local(size);

            // Test Average
            let steps_prim = linkage_with_algo(&m, Method::Average, Algorithm::Primitive);
            let steps_nn = linkage_with_algo(&m, Method::Average, Algorithm::NnChain);

            assert_eq!(steps_prim.len(), size - 1);
            assert_eq!(steps_nn.len(), size - 1);

            for (j, (s1, s2)) in steps_prim.iter().zip(steps_nn.iter()).enumerate() {
                assert!(
                    (s1.distance - s2.distance).abs() < 1e-5,
                    "Iter {}, Size {}, Step {}: Average distance mismatch {} vs {}",
                    i,
                    size,
                    j,
                    s1.distance,
                    s2.distance
                );
            }

            // Test Ward (which uses squared optimization)
            let steps_prim_ward = linkage_with_algo(&m, Method::Ward, Algorithm::Primitive);
            let steps_nn_ward = linkage_with_algo(&m, Method::Ward, Algorithm::NnChain);

            for (j, (s1, s2)) in steps_prim_ward.iter().zip(steps_nn_ward.iter()).enumerate() {
                assert!(
                    (s1.distance - s2.distance).abs() < 1e-5,
                    "Iter {}, Size {}, Step {}: Ward distance mismatch {} vs {}",
                    i,
                    size,
                    j,
                    s1.distance,
                    s2.distance
                );
            }
        }
    }

    #[test]
    fn test_monotonicity() {
        // Check if merge heights are non-decreasing for monotonic methods
        let m = create_random_matrix_local(30);

        let methods = [
            Method::Single,
            Method::Complete,
            Method::Average,
            Method::Weighted,
            Method::Ward,
        ];

        for method in methods {
            let steps = linkage(&m, method);
            for i in 0..steps.len() - 1 {
                assert!(
                    steps[i].distance <= steps[i + 1].distance + 1e-6,
                    "Method {:?} is not monotonic at step {}: {} > {}",
                    method,
                    i,
                    steps[i].distance,
                    steps[i + 1].distance
                );
            }
        }
    }

    #[test]
    fn test_edge_cases() {
        // N=0
        let m0 = NamedMatrix::new(vec![]);
        assert!(linkage(&m0, Method::Average).is_empty());

        // N=1
        let m1 = NamedMatrix::new(vec!["A".to_string()]);
        assert!(linkage(&m1, Method::Average).is_empty());

        // N=2
        let mut m2 = NamedMatrix::new(vec!["A".to_string(), "B".to_string()]);
        m2.set(0, 1, 0.5);
        let steps = linkage(&m2, Method::Average);
        assert_eq!(steps.len(), 1);
        assert_eq!(steps[0].distance, 0.5);
    }

    #[test]
    fn test_to_tree() {
        // Steps from Single Linkage:
        // 1. Merge 0-1 (d=1.0).
        // 2. Merge 3-2 (d=2.0).
        let steps = vec![
            Step {
                cluster1: 0,
                cluster2: 1,
                distance: 1.0,
                size: 2,
            },
            Step {
                cluster1: 3,
                cluster2: 2,
                distance: 2.0,
                size: 3,
            },
        ];
        let names = vec!["A".to_string(), "B".to_string(), "C".to_string()];

        let tree = to_tree(&steps, &names);

        // Root is the last created node
        assert!(tree.get_root().is_some());
        assert_eq!(tree.len(), 5); // 3 leaves + 2 internal

        // Check topology
        // Root should have children 3 and 2
        // Node 3 should have children 0 and 1

        // Find leaf nodes by name
        // Use traversal to find nodes
        let leaf_a_id = tree
            .get_leaves()
            .iter()
            .find(|&&id| tree.get_node(id).unwrap().name.as_deref() == Some("A"))
            .copied()
            .unwrap();
        let leaf_b_id = tree
            .get_leaves()
            .iter()
            .find(|&&id| tree.get_node(id).unwrap().name.as_deref() == Some("B"))
            .copied()
            .unwrap();
        let leaf_c_id = tree
            .get_leaves()
            .iter()
            .find(|&&id| tree.get_node(id).unwrap().name.as_deref() == Some("C"))
            .copied()
            .unwrap();

        let leaf_a = tree.get_node(leaf_a_id).unwrap();
        let leaf_b = tree.get_node(leaf_b_id).unwrap();
        let leaf_c = tree.get_node(leaf_c_id).unwrap();

        // Check heights/lengths
        // Leaf height = 0
        // Node 3 height = 1.0 / 2 = 0.5
        // Root height = 2.0 / 2 = 1.0

        // Length A -> 3: 0.5 - 0 = 0.5
        assert_eq!(leaf_a.length, Some(0.5));
        assert_eq!(leaf_b.length, Some(0.5));

        // Length C -> Root: 1.0 - 0 = 1.0
        assert_eq!(leaf_c.length, Some(1.0));

        // Node 3 -> Root: 1.0 - 0.5 = 0.5
        // We need to find Node 3 (parent of A and B)
        let parent_a = tree.get_node(leaf_a.parent.unwrap()).unwrap();
        assert_eq!(parent_a.length, Some(0.5));
    }
}
