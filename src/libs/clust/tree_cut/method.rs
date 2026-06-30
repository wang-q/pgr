/// Cut strategies for tree partitioning.
#[derive(Debug, Clone, Copy)]
pub enum Method {
    /// Cut into exactly K clusters.
    ///
    /// Iteratively splits the cluster with the largest height (distance to farthest leaf).
    K(usize),

    /// Cut at specific height (distance from leaves).
    ///
    /// Useful for ultrametric trees where height represents time/divergence.
    Height(f64),

    /// Cut at specific distance from root.
    ///
    /// Useful for defining clusters based on divergence from a common ancestor (root).
    RootDist(f64),

    /// TreeCluster: Max pairwise distance in clade <= threshold.
    ///
    /// Ensures that for every cluster, the maximum distance between any two leaves
    /// in that cluster is at most `threshold`.
    MaxClade(f64),

    /// TreeCluster: Average pairwise distance in clade <= threshold.
    AvgClade(f64),

    /// TreeCluster: Median pairwise distance in clade <= threshold.
    MedClade(f64),

    /// TreeCluster: Sum of branch lengths in clade <= threshold.
    SumBranch(f64),

    /// SciPy: Inconsistent coefficient <= threshold.
    ///
    /// Splits nodes if their inconsistency coefficient > threshold.
    /// Requires checking inconsistency of all descendants.
    /// Parameters: (threshold, depth).
    Inconsistent(f64, usize),

    /// TreeCluster: Single Linkage.
    ///
    /// Removes any edge (branch) with length > threshold.
    /// The resulting connected components (subtrees) form clusters.
    /// Note: This is equivalent to `Height` on ultrametric trees but generalizes to any tree.
    /// It effectively breaks "long branches".
    SingleLinkage(f64),
}

/// Supported cut method names, in detection priority order.
/// Excludes `dynamic-tree` and `dynamic-hybrid` which are handled separately.
pub const METHOD_NAMES: &[&str] = &[
    "k",
    "height",
    "root-dist",
    "max-clade",
    "avg-clade",
    "med-clade",
    "sum-branch",
    "leaf-dist-max",
    "leaf-dist-min",
    "leaf-dist-avg",
    "max-edge",
    "inconsistent",
];

/// Build a Method from a name and threshold value.
///
/// For "leaf-dist-*" methods, `leaf_depths` must be provided as `(min, max, avg)`.
pub fn build_method(
    name: &str,
    val: f64,
    deep: usize,
    leaf_depths: Option<(f64, f64, f64)>,
) -> Result<Method, String> {
    match name {
        "k" => Ok(Method::K(val as usize)),
        "height" => Ok(Method::Height(val)),
        "root-dist" => Ok(Method::RootDist(val)),
        "max-clade" => Ok(Method::MaxClade(val)),
        "avg-clade" => Ok(Method::AvgClade(val)),
        "med-clade" => Ok(Method::MedClade(val)),
        "sum-branch" => Ok(Method::SumBranch(val)),
        "leaf-dist-max" => leaf_depths
            .map(|(_, max, _)| Method::RootDist(max - val))
            .ok_or_else(|| "leaf depths required for leaf-dist-max".to_string()),
        "leaf-dist-min" => leaf_depths
            .map(|(min, _, _)| Method::RootDist(min - val))
            .ok_or_else(|| "leaf depths required for leaf-dist-min".to_string()),
        "leaf-dist-avg" => leaf_depths
            .map(|(_, _, avg)| Method::RootDist(avg - val))
            .ok_or_else(|| "leaf depths required for leaf-dist-avg".to_string()),
        "max-edge" => Ok(Method::SingleLinkage(val)),
        "inconsistent" => Ok(Method::Inconsistent(val, deep)),
        _ => Err(format!("unknown method: {}", name)),
    }
}
