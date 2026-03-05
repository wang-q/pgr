use super::Partition;
use crate::libs::pairmat::NamedMatrix;
use crate::libs::phylo::tree::stat::compute_node_heights;
use crate::libs::phylo::tree::traversal::postorder;
use crate::libs::phylo::tree::Tree;
use std::collections::HashMap;

pub struct HybridOptions {
    pub min_cluster_size: usize,
    pub dist_matrix: NamedMatrix,
    pub cut_height: Option<f64>,
    pub deep_split: usize,             // 0..4, default 1
    pub max_core_scatter: Option<f64>, // relative 0..1
    pub min_gap: Option<f64>,          // relative 0..1
    pub pam_stage: bool,
    pub pam_respects_dendro: bool,
    pub max_pam_dist: Option<f64>,
    pub respect_small_clusters: bool,
}

struct BranchInfo {
    members: Vec<usize>, // Matrix indices
    is_basic: bool,
    scatter: f64,
}

pub fn cutree_hybrid(tree: &Tree, options: HybridOptions) -> anyhow::Result<Partition> {
    let matrix = &options.dist_matrix;
    let min_size = options.min_cluster_size;

    // 1. Calculate Cut Parameters (CutHeight, RefHeight)
    let node_heights = compute_node_heights(tree);

    // Get all merge heights (internal node heights)
    let mut merge_heights: Vec<f64> = node_heights
        .iter()
        .filter(|(&id, _)| {
            if let Some(node) = tree.get_node(id) {
                !node.is_leaf()
            } else {
                false
            }
        })
        .map(|(_, &h)| h)
        .collect();
    merge_heights.sort_by(|a, b| a.partial_cmp(b).unwrap());

    if merge_heights.is_empty() {
        // Single node tree?
        return Ok(Partition {
            num_clusters: 0,
            assignment: HashMap::new(),
        });
    }

    let n_merge = merge_heights.len();
    let ref_idx = (n_merge as f64 * 0.05).round() as usize;
    let ref_idx = ref_idx.max(1).min(n_merge - 1); // 5th percentile
    let ref_height = merge_heights[ref_idx];
    let max_height = merge_heights.last().copied().unwrap_or(0.0);

    let cut_height = options
        .cut_height
        .unwrap_or_else(|| 0.99 * (max_height - ref_height) + ref_height);

    // 2. Calculate Absolute Thresholds
    let deep_split = options.deep_split.clamp(0, 4);
    let def_mcs = [0.64, 0.73, 0.82, 0.91, 0.95];
    let def_mg: Vec<f64> = def_mcs.iter().map(|&x| (1.0 - x) * 0.75).collect();

    let rel_max_scatter = options.max_core_scatter.unwrap_or(def_mcs[deep_split]);
    let rel_min_gap = options.min_gap.unwrap_or(def_mg[deep_split]);

    let abs_max_scatter = ref_height + rel_max_scatter * (cut_height - ref_height);
    let abs_min_gap = rel_min_gap * (cut_height - ref_height);

    // 3. Map Leaves to Matrix Indices
    let mut node_to_mat_idx = HashMap::new();
    let mut mat_idx_to_node = HashMap::new();
    for leaf_id in tree.get_leaves() {
        if let Some(node) = tree.get_node(leaf_id) {
            if let Some(name) = &node.name {
                if let Some(idx) = matrix.get_index(name) {
                    node_to_mat_idx.insert(leaf_id, idx);
                    mat_idx_to_node.insert(idx, leaf_id);
                }
            }
        }
    }

    // 4. Bottom-Up Traversal (Core Detection)
    let mut branch_infos: HashMap<usize, BranchInfo> = HashMap::new();
    let mut final_clusters: Vec<Vec<usize>> = Vec::new();
    let mut small_clusters: Vec<Vec<usize>> = Vec::new(); // Clusters that failed size check but passed shape

    let root_id = match tree.get_root() {
        Some(r) => r,
        None => {
            return Ok(Partition {
                num_clusters: 0,
                assignment: HashMap::new(),
            });
        }
    };

    let post_order = postorder(tree, root_id);

    for &node_id in &post_order {
        let node = tree.get_node(node_id).unwrap();

        if node.is_leaf() {
            if let Some(&mat_idx) = node_to_mat_idx.get(&node_id) {
                branch_infos.insert(
                    node_id,
                    BranchInfo {
                        members: vec![mat_idx],
                        is_basic: true,
                        scatter: 0.0,
                    },
                );
            }
            continue;
        }

        // Internal Node
        let mut child_infos = Vec::new();
        let mut all_members = Vec::new();
        for &child_id in &node.children {
            if let Some(info) = branch_infos.remove(&child_id) {
                all_members.extend_from_slice(&info.members);
                child_infos.push(info);
            }
        }

        if child_infos.is_empty() {
            continue;
        }

        let node_height = *node_heights.get(&node_id).unwrap_or(&0.0);

        let current_scatter = if all_members.len() >= min_size {
            calculate_core_scatter(&all_members, matrix, min_size)
        } else {
            0.0
        };

        // Determine Merge Logic
        let do_merge;

        if node_height > cut_height {
            do_merge = false;
        } else {
            let mut all_basic = true;
            let mut any_fail = false;

            for info in &child_infos {
                if !info.is_basic {
                    all_basic = false;
                } else {
                    let gap = node_height - info.scatter;
                    let fail_shape = info.scatter > abs_max_scatter || gap < abs_min_gap;
                    let fail_size = info.members.len() < min_size;

                    if fail_shape || fail_size {
                        any_fail = true;
                    }
                }
            }

            if all_basic {
                // If all children are basic, we merge if ANY of them fails (size or shape).
                // This logic tries to grow clusters up until they break criteria.
                if any_fail {
                    do_merge = true;
                } else {
                    // If all children satisfy criteria, we do NOT merge, effectively stopping
                    // at the smallest valid clusters (matches R's bottom-up behavior).
                    do_merge = false;
                }
            } else {
                // Mixed Basic and Composite.
                do_merge = false;
            }
        }

        if do_merge {
            // Form Basic Cluster
            branch_infos.insert(
                node_id,
                BranchInfo {
                    members: all_members,
                    is_basic: true,
                    scatter: current_scatter,
                },
            );
        } else {
            // Form Composite Cluster
            // Finalize Strong Basic children
            for info in child_infos {
                if info.is_basic {
                    let gap = node_height - info.scatter;
                    let fail_shape = info.scatter > abs_max_scatter || gap < abs_min_gap;
                    let fail_size = info.members.len() < min_size;

                    if !fail_shape {
                        if !fail_size {
                            final_clusters.push(info.members);
                        } else if options.respect_small_clusters {
                            small_clusters.push(info.members);
                        }
                    }
                }
                // Composite children are already handled (their sub-clusters finalized).
            }

            branch_infos.insert(
                node_id,
                BranchInfo {
                    members: all_members,
                    is_basic: false,
                    scatter: 0.0,
                },
            );
        }
    }

    // Handle Root (or remaining top node)
    if let Some(root_id) = tree.get_root() {
        if let Some(info) = branch_infos.remove(&root_id) {
            if info.is_basic {
                let fail_shape = info.scatter > abs_max_scatter; // Gap is undefined at root/top? Or check against cut_height?
                let fail_size = info.members.len() < min_size;

                if !fail_shape {
                    if !fail_size {
                        final_clusters.push(info.members);
                    } else if options.respect_small_clusters {
                        small_clusters.push(info.members);
                    }
                }
            }
        }
    }

    // 5. Initial Assignment & Medoid Calculation
    let mut assignment = HashMap::new();
    let mut medoids: HashMap<usize, usize> = HashMap::new(); // ClusterID -> MedoidIdx
    let mut medoid_node_ids: HashMap<usize, usize> = HashMap::new(); // ClusterID -> NodeId (for LCA check)

    // Assign initial clusters
    for (i, members) in final_clusters.iter().enumerate() {
        let cid = i + 1;

        // Find Medoid
        let mut min_total_dist = f64::MAX;
        let mut best_medoid = members[0];

        for &candidate in members {
            let mut total_dist = 0.0;
            for &other in members {
                if candidate != other {
                    total_dist += matrix.get(candidate, other) as f64;
                }
            }
            if total_dist < min_total_dist {
                min_total_dist = total_dist;
                best_medoid = candidate;
            }
        }
        medoids.insert(cid, best_medoid);
        if let Some(&nid) = mat_idx_to_node.get(&best_medoid) {
            medoid_node_ids.insert(cid, nid);
        }

        for &idx in members {
            if let Some(&node_id) = mat_idx_to_node.get(&idx) {
                assignment.insert(node_id, cid);
            }
        }
    }

    // 6. PAM Stage
    if options.pam_stage && !medoids.is_empty() {
        let max_pam_dist = options.max_pam_dist.unwrap_or(cut_height);

        // Identify unassigned nodes (Cluster 0).
        // If respect_small_clusters is true, we first try to assign small clusters as blocks
        // before handling individual objects.

        // If respect_small_clusters is true, we first try to assign small_clusters as blocks
        if options.respect_small_clusters {
            for members in &small_clusters {
                // Find best cluster for this block
                let mut best_cid = 0;
                let mut min_avg_dist = f64::MAX;

                for (&cid, &medoid_idx) in &medoids {
                    let mut total_dist = 0.0;
                    for &m in members {
                        total_dist += matrix.get(m, medoid_idx) as f64;
                    }
                    let avg_dist = total_dist / members.len() as f64;

                    if avg_dist < min_avg_dist {
                        // Check dendro constraint
                        let mut valid = true;
                        if options.pam_respects_dendro {
                            if let Some(&medoid_node) = medoid_node_ids.get(&cid) {
                                // Check LCA for all members? Or just one?
                                // R says "an object". For a cluster, maybe check representative?
                                // Let's check the first member (approximation) or all.
                                // Strict: All members must be compatible.
                                for &m in members {
                                    if let Some(&m_node) = mat_idx_to_node.get(&m) {
                                        if let Ok(lca) =
                                            tree.get_common_ancestor(&m_node, &medoid_node)
                                        {
                                            let h = *node_heights.get(&lca).unwrap_or(&0.0);
                                            if h > cut_height {
                                                valid = false;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        if valid {
                            min_avg_dist = avg_dist;
                            best_cid = cid;
                        }
                    }
                }

                if min_avg_dist <= max_pam_dist {
                    for &idx in members {
                        if let Some(&node_id) = mat_idx_to_node.get(&idx) {
                            assignment.insert(node_id, best_cid);
                        }
                    }
                }
            }
        }

        // Now assign remaining individual objects (singleton outliers).
        // If respect_small_clusters is TRUE, we skip objects already assigned as part of a block.
        // This ensures we don't break the structure of small clusters identified in the first stage.

        let all_indices: Vec<usize> = node_to_mat_idx.values().cloned().collect();
        for &idx in &all_indices {
            let node_id = *mat_idx_to_node.get(&idx).unwrap();
            if assignment.contains_key(&node_id) {
                continue;
            }

            let mut best_cid = 0;
            let mut min_dist = f64::MAX;

            for (&cid, &medoid_idx) in &medoids {
                let d = matrix.get(idx, medoid_idx) as f64;

                if d < min_dist {
                    // Check dendro constraint
                    let mut valid = true;
                    if options.pam_respects_dendro {
                        if let Some(&medoid_node) = medoid_node_ids.get(&cid) {
                            if let Ok(lca) = tree.get_common_ancestor(&node_id, &medoid_node) {
                                let h = *node_heights.get(&lca).unwrap_or(&0.0);
                                if h > cut_height {
                                    valid = false;
                                }
                            }
                        }
                    }

                    if valid {
                        min_dist = d;
                        best_cid = cid;
                    }
                }
            }

            if min_dist <= max_pam_dist {
                assignment.insert(node_id, best_cid);
            }
        }
    }

    // Fill missing nodes with 0
    for node_id in tree.get_leaves() {
        assignment.entry(node_id).or_insert(0);
    }

    let _num_clusters = medoids.keys().len(); // Approximate, some might be empty now?
                                              // Recalculate strict num_clusters
    let max_cid = assignment
        .values()
        .filter(|&&v| v > 0)
        .max()
        .copied()
        .unwrap_or(0);

    Ok(Partition {
        num_clusters: max_cid,
        assignment,
    })
}

fn calculate_core_scatter(members: &[usize], matrix: &NamedMatrix, min_cluster_size: usize) -> f64 {
    let n = members.len();
    if n == 0 {
        return 0.0;
    }

    // Calculate average distance for each point to others
    let mut point_avg_dists = Vec::with_capacity(n);
    for &i in members {
        let mut sum = 0.0;
        for &j in members {
            if i != j {
                sum += matrix.get(i, j) as f64;
            }
        }
        point_avg_dists.push((sum / (n - 1).max(1) as f64, i));
    }

    point_avg_dists.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let base_core_size = (min_cluster_size as f64 / 2.0 + 1.0) as usize;
    let eff_core_size = if base_core_size < n {
        (base_core_size as f64 + (n - base_core_size) as f64).sqrt() as usize + base_core_size
    } else {
        n
    };

    let core_members: Vec<usize> = point_avg_dists
        .iter()
        .take(eff_core_size)
        .map(|x| x.1)
        .collect();

    // Calculate average distance within core
    let mut core_sum = 0.0;
    let mut count = 0;
    for i in 0..core_members.len() {
        for j in (i + 1)..core_members.len() {
            core_sum += matrix.get(core_members[i], core_members[j]) as f64;
            count += 1;
        }
    }

    if count == 0 {
        0.0
    } else {
        core_sum / count as f64
    }
}
