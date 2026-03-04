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

/// Perform hierarchical clustering on a condensed distance matrix.
///
/// Returns a list of steps (merges) forming the dendrogram.
/// The length of the result will be N-1 for a matrix of size N.
pub fn linkage(matrix: &NamedMatrix, method: Method) -> Vec<Step> {
    match method {
        Method::Single => linkage_primitive(matrix, method), // Should use MST in Phase 2
        _ => linkage_primitive(matrix, method),              // Use primitive O(N^3) for Phase 1
    }
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
    let mut cluster_to_node: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();

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
        let h_parent = step.distance;

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
fn linkage_primitive(matrix: &NamedMatrix, method: Method) -> Vec<Step> {
    let n = matrix.size();
    if n < 2 {
        return vec![];
    }

    // Work on a mutable copy of distances
    let mut condensed = CondensedMatrix::from_vec(n, matrix.values().to_vec());

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

        steps.push(Step {
            cluster1: id1,
            cluster2: id2,
            distance: min_dist,
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

            let new_dist = lance_williams(
                method, d_uk, d_vk, d_uv, size1, size2, size[k],
            );
            
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
            // d(u+v, k)^2 = (n_u * d(u,k)^2 + n_v * d(v,k)^2 - n_u*n_v*d(u,v)^2 / n_uv) / n_uv
            // Implemented using squared distances logic for consistency with Ward/Median.
            
            let d2_uk = d_uk * d_uk;
            let d2_vk = d_vk * d_vk;
            let d2_uv = d_uv * d_uv;
            
            let d2_new = (n_u * d2_uk + n_v * d2_vk - (n_u * n_v * d2_uv) / n_uv) / n_uv;
            if d2_new < 0.0 { 0.0 } else { d2_new.sqrt() }
        }
        Method::Median => {
            // WPGMC (Weighted Pair Group Method with Centroid Averaging)
            // d_new = 0.5*d_uk + 0.5*d_vk - 0.25*d_uv (Squared)
            let d2_uk = d_uk * d_uk;
            let d2_vk = d_vk * d_vk;
            let d2_uv = d_uv * d_uv;
            let d2_new = 0.5 * d2_uk + 0.5 * d2_vk - 0.25 * d2_uv;
            if d2_new < 0.0 { 0.0 } else { d2_new.sqrt() }
        }
        Method::Ward => {
            // Ward's method (minimal variance) - specifically Ward.D2
            // d_new^2 = ((n_u+n_k)*d_uk^2 + (n_v+n_k)*d_vk^2 - n_k*d_uv^2) / n_uvk
            let d2_uk = d_uk * d_uk;
            let d2_vk = d_vk * d_vk;
            let d2_uv = d_uv * d_uv;
            
            let d2_new = ((n_u + n_k) * d2_uk + (n_v + n_k) * d2_vk - n_k * d2_uv) / n_uvk;
            if d2_new < 0.0 { 0.0 } else { d2_new.sqrt() }
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

        assert_eq!(steps[1].cluster1, 3); // 3 is the new id for {0,1}
        assert_eq!(steps[1].cluster2, 2);
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
    fn test_to_tree() {
        // Steps from Single Linkage:
        // 1. Merge 0-1 (d=1.0).
        // 2. Merge 3-2 (d=2.0).
        let steps = vec![
            Step { cluster1: 0, cluster2: 1, distance: 1.0, size: 2 },
            Step { cluster1: 3, cluster2: 2, distance: 2.0, size: 3 },
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
        let leaf_a_id = tree.get_leaves().iter().find(|&&id| {
            tree.get_node(id).unwrap().name.as_deref() == Some("A")
        }).copied().unwrap();
        let leaf_b_id = tree.get_leaves().iter().find(|&&id| {
            tree.get_node(id).unwrap().name.as_deref() == Some("B")
        }).copied().unwrap();
        let leaf_c_id = tree.get_leaves().iter().find(|&&id| {
            tree.get_node(id).unwrap().name.as_deref() == Some("C")
        }).copied().unwrap();

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
