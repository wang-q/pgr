use crate::libs::phylo::node::NodeId;
use crate::libs::phylo::tree::Tree;
use std::collections::HashMap;

pub mod clade;
pub mod dynamic;
pub mod hybrid;
pub mod inconsistent;
pub mod link;
pub mod method;
pub mod simple;

pub use method::Method;

/// Result of a cut operation.
pub struct Partition {
    /// Map from Leaf NodeId to Cluster ID (1-based).
    pub assignment: HashMap<NodeId, usize>,
    /// Total number of clusters formed.
    pub num_clusters: usize,
}

impl Partition {
    /// Create a new empty partition.
    pub fn new() -> Self {
        Self {
            assignment: HashMap::new(),
            num_clusters: 0,
        }
    }

    /// Get members of each cluster.
    ///
    /// Returns a map where keys are Cluster IDs (1-based) and values are lists of Leaf NodeIds.
    pub fn get_clusters(&self) -> HashMap<usize, Vec<NodeId>> {
        let mut clusters = HashMap::new();
        for (&node_id, &cluster_id) in &self.assignment {
            clusters
                .entry(cluster_id)
                .or_insert_with(Vec::new)
                .push(node_id);
        }
        clusters
    }

    /// Compute summary statistics for the partition.
    /// Returns (num_clusters, num_singletons, num_non_singletons, max_cluster_size).
    pub fn get_stats(&self) -> (usize, usize, usize, usize) {
        let mut sizes = HashMap::new();
        for &cluster_id in self.assignment.values() {
            *sizes.entry(cluster_id).or_insert(0) += 1;
        }
        let mut singletons = 0;
        let mut non_singletons = 0;
        let mut max_size = 0;
        for &size in sizes.values() {
            if size == 1 {
                singletons += 1;
            } else {
                non_singletons += 1;
            }
            if size > max_size {
                max_size = size;
            }
        }
        (self.num_clusters, singletons, non_singletons, max_size)
    }
}

impl Default for Partition {
    fn default() -> Self {
        Self::new()
    }
}

/// Cut the tree according to the specified method.
///
/// # Arguments
///
/// *   `tree` - The phylogenetic tree to cut.
/// *   `method` - The cutting method/strategy.
///
/// # Returns
///
/// A `Result` containing the `Partition` or an error message.
pub fn cut(tree: &Tree, method: Method) -> Result<Partition, String> {
    if tree.is_empty() {
        return Ok(Partition::new());
    }

    match method {
        Method::K(k) => simple::cut_k(tree, k),
        Method::Height(h) => simple::cut_height(tree, h),
        Method::RootDist(d) => simple::cut_root_dist(tree, d),
        Method::MaxClade(t) => clade::cut_max_clade(tree, t),
        Method::AvgClade(t) => clade::cut_avg_clade(tree, t),
        Method::MedClade(t) => clade::cut_med_clade(tree, t),
        Method::SumBranch(t) => clade::cut_sum_branch(tree, t),
        Method::Inconsistent(t, d) => inconsistent::cut_inconsistent(tree, t, d),
        Method::SingleLinkage(t) => link::cut_single_linkage(tree, t),
    }
}

// --- Helpers ---

/// Compute max distance from each node to its leaves
pub(crate) fn compute_heights(tree: &Tree, root: NodeId) -> Result<HashMap<NodeId, f64>, String> {
    let mut heights = HashMap::new();
    let post_order = crate::libs::phylo::tree::traversal::postorder(tree, root);

    for id in post_order {
        let node = tree.get_node(id).unwrap();
        if node.children.is_empty() {
            heights.insert(id, 0.0);
        } else {
            let mut max_h = 0.0;
            for &child in &node.children {
                let child_node = tree.get_node(child).unwrap();
                let len = child_node.length.unwrap_or(0.0); // If None, assume 0
                let h = heights[&child] + len;
                if h > max_h {
                    max_h = h;
                }
            }
            heights.insert(id, max_h);
        }
    }
    Ok(heights)
}

/// Assign leaves to clusters based on cluster roots
pub(crate) fn assign_clusters(
    tree: &Tree,
    cluster_roots: Vec<NodeId>,
) -> Result<Partition, String> {
    let mut part = Partition::new();
    let mut cluster_id = 0;

    for root in cluster_roots {
        cluster_id += 1;
        let leaves = crate::libs::phylo::tree::stat::get_leaves(tree, root);
        for leaf in leaves {
            part.assignment.insert(leaf, cluster_id);
        }
    }

    part.num_clusters = cluster_id;
    Ok(part)
}

/// Mask low-support internal nodes by setting branch length to infinity (TreeCluster semantics).
pub fn apply_support_filter(tree: &mut Tree, threshold: f64) {
    let len = tree.len();
    for i in 0..len {
        let should_mask = {
            if let Some(node) = tree.get_node(i) {
                if !node.children.is_empty() {
                    let support = node
                        .name
                        .as_ref()
                        .and_then(|n| n.parse::<f64>().ok())
                        .unwrap_or(100.0);
                    support < threshold
                } else {
                    false
                }
            } else {
                false
            }
        };

        if should_mask {
            if let Some(node) = tree.get_node_mut(i) {
                node.length = Some(f64::INFINITY);
            }
        }
    }
}

// --- Cluster output helpers ---

/// Representative selection mode for clusters.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RepMode {
    /// Member closest to root (alphabetical tie-break).
    Root,
    /// Alphabetically first member.
    First,
    /// Member with min sum of distances to others (alphabetical tie-break).
    Medoid,
}

impl RepMode {
    /// Parse a rep mode from a string ("root", "first", "medoid").
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "root" => Ok(RepMode::Root),
            "first" => Ok(RepMode::First),
            "medoid" => Ok(RepMode::Medoid),
            _ => Err(format!("unsupported rep method: {}", s)),
        }
    }
}

/// A cluster of tree leaves with members sorted alphabetically by name.
#[derive(Debug, Clone)]
pub struct Cluster {
    /// (NodeId, name) pairs, sorted alphabetically by name.
    pub members: Vec<(NodeId, String)>,
    /// Index of the representative in `members` (None if cluster is empty).
    pub rep_index: Option<usize>,
}

impl Cluster {
    /// Get the representative name, if any.
    pub fn rep_name(&self) -> Option<&str> {
        self.rep_index
            .and_then(|i| self.members.get(i).map(|(_, n)| n.as_str()))
    }
}

/// Select the representative index for a cluster.
/// Returns `Some(index)` or `None` if the cluster is empty.
pub fn find_representative(
    cluster: &Cluster,
    tree: &Tree,
    rep_mode: RepMode,
    root_dists: &HashMap<NodeId, f64>,
) -> Option<usize> {
    let members = &cluster.members;
    if members.is_empty() {
        return None;
    }
    match rep_mode {
        RepMode::First => Some(0),
        RepMode::Root => members
            .iter()
            .enumerate()
            .min_by(|(_, (id1, _)), (_, (id2, _))| {
                let d1 = root_dists.get(id1).copied().unwrap_or(f64::MAX);
                let d2 = root_dists.get(id2).copied().unwrap_or(f64::MAX);
                d1.partial_cmp(&d2).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(idx, _)| idx),
        RepMode::Medoid => {
            let ids: Vec<NodeId> = members.iter().map(|(id, _)| *id).collect();
            crate::libs::phylo::tree::query::tree_medoid(tree, &ids)
        }
    }
}

/// Convert a partition into clusters with representatives selected.
/// Clusters are sorted by size (descending), then by first member name.
pub fn partition_to_clusters(
    partition: &Partition,
    tree: &Tree,
    rep_mode: RepMode,
) -> Vec<Cluster> {
    let clusters_map = partition.get_clusters();
    let root_dists = crate::libs::phylo::tree::stat::compute_root_distances(tree);

    let mut clusters: Vec<Cluster> = Vec::with_capacity(clusters_map.len());
    for members in clusters_map.values() {
        let mut member_info: Vec<(NodeId, String)> = Vec::with_capacity(members.len());
        for &mid in members {
            if let Some(node) = tree.get_node(mid) {
                let name = node.name.clone().unwrap_or_else(|| format!("Leaf_{}", mid));
                member_info.push((mid, name));
            }
        }
        member_info.sort_by(|a, b| a.1.cmp(&b.1));

        let mut cluster = Cluster {
            members: member_info,
            rep_index: None,
        };
        cluster.rep_index = find_representative(&cluster, tree, rep_mode, &root_dists);
        clusters.push(cluster);
    }

    // Sort clusters: first by size (descending), then by first member name.
    clusters.sort_by(|a, b| match b.members.len().cmp(&a.members.len()) {
        std::cmp::Ordering::Equal => {
            let name_a = a.members.first().map(|s| s.1.as_str()).unwrap_or("");
            let name_b = b.members.first().map(|s| s.1.as_str()).unwrap_or("");
            name_a.cmp(name_b)
        }
        other => other,
    });

    clusters
}

/// Format clusters into output string.
/// `format` must be "cluster" or "pair".
pub fn format_clusters(clusters: &[Cluster], format: &str) -> Result<String, String> {
    let mut out = String::new();
    match format {
        "cluster" => {
            for c in clusters {
                if let Some(rep_idx) = c.rep_index {
                    let mut names: Vec<&str> = c.members.iter().map(|(_, n)| n.as_str()).collect();
                    if rep_idx != 0 {
                        names.swap(0, rep_idx);
                        names[1..].sort();
                    }
                    out.push_str(&names.join("\t"));
                    out.push('\n');
                }
            }
        }
        "pair" => {
            for c in clusters {
                if let Some(rep_name) = c.rep_name() {
                    for (_, member_name) in &c.members {
                        out.push_str(rep_name);
                        out.push('\t');
                        out.push_str(member_name);
                        out.push('\n');
                    }
                }
            }
        }
        _ => return Err(format!("unsupported output format: {}", format)),
    }
    Ok(out)
}

/// Format a partition as scan-mode TSV rows.
/// Each row is "group_label\tcluster_id\tmember_name", clusters ordered by ID.
pub fn format_scan_rows(partition: &Partition, tree: &Tree, group_label: &str) -> String {
    let clusters_map = partition.get_clusters();
    let mut cluster_ids: Vec<usize> = clusters_map.keys().copied().collect();
    cluster_ids.sort();

    let mut out = String::new();
    for (i, cid) in cluster_ids.iter().enumerate() {
        let cluster_label = i + 1;
        if let Some(members) = clusters_map.get(cid) {
            let mut member_names: Vec<String> = Vec::with_capacity(members.len());
            for &mid in members {
                if let Some(node) = tree.get_node(mid) {
                    let name = node.name.clone().unwrap_or_else(|| format!("Leaf_{}", mid));
                    member_names.push(name);
                }
            }
            member_names.sort();
            for name in member_names {
                out.push_str(&format!("{}\t{}\t{}\n", group_label, cluster_label, name));
            }
        }
    }
    out
}
