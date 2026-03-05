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
