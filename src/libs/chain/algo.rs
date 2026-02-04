// use std::cmp::Ordering;

/// Trait representing an item that can be indexed in the KD-tree for chaining.
///
/// Each item corresponds to an alignment block with query and target coordinates
/// and a score.
pub trait ChainItem {
    /// Query sequence start coordinate (0-based, inclusive).
    fn q_start(&self) -> u64;
    /// Query sequence end coordinate (0-based, exclusive).
    fn q_end(&self) -> u64;
    /// Target sequence start coordinate (0-based, inclusive).
    fn t_start(&self) -> u64;
    /// Target sequence end coordinate (0-based, exclusive).
    fn t_end(&self) -> u64;
    /// Score of the alignment block.
    fn score(&self) -> f64;
}

/// A leaf node in the KD-tree storing a reference to an item.
pub struct KdLeaf<T> {
    pub item: T,
    pub original_idx: usize,
    pub best_pred_idx: Option<usize>,
    pub total_score: f64,
    pub hit: bool,
}

/// A node in the KD-tree, either a Leaf or an Internal node.
pub enum KdNode {
    Leaf {
        leaf_idx: usize,
        max_q: u64,
        max_t: u64,
        max_score: f64,
    },
    Internal {
        cut_coord: u64,
        lo: Box<KdNode>,
        hi: Box<KdNode>,
        max_q: u64,
        max_t: u64,
        max_score: f64,
    },
}

impl KdNode {
    pub fn max_q(&self) -> u64 {
        match self {
            KdNode::Leaf { max_q, .. } => *max_q,
            KdNode::Internal { max_q, .. } => *max_q,
        }
    }
    pub fn max_t(&self) -> u64 {
        match self {
            KdNode::Leaf { max_t, .. } => *max_t,
            KdNode::Internal { max_t, .. } => *max_t,
        }
    }
    pub fn max_score(&self) -> f64 {
        match self {
            KdNode::Leaf { max_score, .. } => *max_score,
            KdNode::Internal { max_score, .. } => *max_score,
        }
    }
}

/// A 2D KD-tree for efficient range queries and predecessor search in chaining algorithms.
///
/// It indexes items based on their query (`q_start`) and target (`t_start`) coordinates.
pub struct KdTree {
    root: Option<Box<KdNode>>,
}

impl KdTree {
    /// Builds a KD-tree from a list of indices into the `items` slice.
    ///
    /// The tree construction alternates between splitting on query and target coordinates.
    pub fn build<T: ChainItem>(indices: &mut [usize], items: &[T]) -> Self {
        if indices.is_empty() {
            return KdTree { root: None };
        }
        KdTree {
            root: Some(Self::build_recursive(indices, items, 0)),
        }
    }

    fn build_recursive<T: ChainItem>(indices: &mut [usize], items: &[T], dim: usize) -> Box<KdNode> {
        if indices.len() == 1 {
            let idx = indices[0];
            let item = &items[idx];
            return Box::new(KdNode::Leaf {
                leaf_idx: idx,
                max_q: item.q_end(),
                max_t: item.t_end(),
                max_score: 0.0,
            });
        }

        if dim == 0 {
            indices.sort_by_key(|&i| items[i].q_start());
        } else {
            indices.sort_by_key(|&i| items[i].t_start());
        }

        let mid = indices.len() / 2;
        let cut_coord = if dim == 0 {
            items[indices[mid]].q_start()
        } else {
            items[indices[mid]].t_start()
        };

        let (left_indices, right_indices) = indices.split_at_mut(mid);

        let lo = Self::build_recursive(left_indices, items, 1 - dim);
        let hi = Self::build_recursive(right_indices, items, 1 - dim);

        let max_q = std::cmp::max(lo.max_q(), hi.max_q());
        let max_t = std::cmp::max(lo.max_t(), hi.max_t());

        Box::new(KdNode::Internal {
            cut_coord,
            lo,
            hi,
            max_q,
            max_t,
            max_score: 0.0,
        })
    }

    /// Updates the max score of a leaf node and propagates the change up the tree.
    ///
    /// This is called during dynamic programming when a better score is found for a chain ending at `leaf_idx`.
    pub fn update_scores<T: ChainItem>(&mut self, leaf_idx: usize, score: f64, items: &[T]) {
        if let Some(root) = &mut self.root {
            Self::update_recursive(root, leaf_idx, score, items, 0);
        }
    }

    fn update_recursive<T: ChainItem>(
        node: &mut KdNode,
        target_idx: usize,
        score: f64,
        items: &[T],
        dim: usize,
    ) {
        match node {
            KdNode::Leaf {
                leaf_idx,
                max_score,
                ..
            } => {
                if *leaf_idx == target_idx {
                    if score > *max_score {
                        *max_score = score;
                    }
                }
            }
            KdNode::Internal {
                cut_coord,
                lo,
                hi,
                max_score,
                ..
            } => {
                if score > *max_score {
                    *max_score = score;
                }

                let coord = if dim == 0 {
                    items[target_idx].q_start()
                } else {
                    items[target_idx].t_start()
                };

                if coord < *cut_coord {
                    Self::update_recursive(lo, target_idx, score, items, 1 - dim);
                } else {
                    Self::update_recursive(hi, target_idx, score, items, 1 - dim);
                }
            }
        }
    }

    /// Finds the best predecessor for the item at `target_idx`.
    ///
    /// It searches the KD-tree for a "chainable" item that maximizes the total score.
    ///
    /// # Arguments
    ///
    /// * `target_idx` - Index of the current item in the `items` slice.
    /// * `current_score` - Current best score for the target item (e.g., just its own score).
    /// * `items` - Slice of all chain items.
    /// * `cost_func` - A closure `Fn(candidate_idx, target_idx) -> Option<new_total_score>`.
    ///   It calculates the score if `candidate` precedes `target`. Returns `None` if they cannot be chained.
    /// * `lower_bound_func` - A closure `Fn(dq, dt) -> lower_bound_cost`.
    ///   It returns a lower bound on the gap cost based on distance in query (`dq`) and target (`dt`) coordinates.
    ///   Used for pruning the search.
    ///
    /// # Returns
    ///
    /// A tuple `(best_score, best_pred_idx)`, where `best_pred_idx` is `Some(index)` if a predecessor was found.
    pub fn best_predecessor<T, F, L>(
        &self,
        target_idx: usize,
        current_score: f64,
        items: &[T],
        cost_func: &F,
        lower_bound_func: &L,
    ) -> (f64, Option<usize>)
    where
        T: ChainItem,
        F: Fn(usize, usize) -> Option<f64>, // (candidate_idx, target_idx) -> Option<new_total_score>
        L: Fn(u64, u64) -> f64, // (dq, dt) -> lower_bound_cost
    {
        let mut best_score = current_score;
        let mut best_pred = None;

        if let Some(root) = &self.root {
            let res = Self::best_recursive(
                root,
                target_idx,
                items,
                cost_func,
                lower_bound_func,
                0,
                best_score,
                best_pred,
            );
            best_score = res.0;
            best_pred = res.1;
        }
        (best_score, best_pred)
    }

    fn best_recursive<T, F, L>(
        node: &KdNode,
        target_idx: usize,
        items: &[T],
        cost_func: &F,
        lower_bound_func: &L,
        dim: usize,
        mut best_score: f64,
        mut best_pred: Option<usize>,
    ) -> (f64, Option<usize>)
    where
        T: ChainItem,
        F: Fn(usize, usize) -> Option<f64>,
        L: Fn(u64, u64) -> f64,
    {
        let target_item = &items[target_idx];
        let node_max_score = node.max_score();
        
        // Pruning 1: Even with 0 cost, can we beat best_score?
        // score = candidate_total + target_score - cost
        // max_possible = node_max + target_score
        if node_max_score + target_item.score() < best_score {
            return (best_score, best_pred);
        }

        // Pruning 2: Geometric distance check
        let dq = if target_item.q_start() > node.max_q() {
            target_item.q_start() - node.max_q()
        } else {
            0
        };
        let dt = if target_item.t_start() > node.max_t() {
            target_item.t_start() - node.max_t()
        } else {
            0
        };
        
        let cost = lower_bound_func(dq, dt);
        if node_max_score + target_item.score() - cost < best_score {
            return (best_score, best_pred);
        }

        match node {
            KdNode::Leaf { leaf_idx, .. } => {
                let cand_idx = *leaf_idx;
                if let Some(new_score) = cost_func(cand_idx, target_idx) {
                    if new_score > best_score {
                        best_score = new_score;
                        best_pred = Some(cand_idx);
                    }
                }
            }
            KdNode::Internal { cut_coord, lo, hi, .. } => {
                let dim_coord = if dim == 0 {
                    target_item.q_start()
                } else {
                    target_item.t_start()
                };

                if dim_coord > *cut_coord {
                    let res = Self::best_recursive(hi, target_idx, items, cost_func, lower_bound_func, 1 - dim, best_score, best_pred);
                    best_score = res.0;
                    best_pred = res.1;
                    
                    let res = Self::best_recursive(lo, target_idx, items, cost_func, lower_bound_func, 1 - dim, best_score, best_pred);
                    best_score = res.0;
                    best_pred = res.1;
                } else {
                    let res = Self::best_recursive(lo, target_idx, items, cost_func, lower_bound_func, 1 - dim, best_score, best_pred);
                    best_score = res.0;
                    best_pred = res.1;
                }
            }
        }
        (best_score, best_pred)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct TestItem {
        q_start: u64,
        q_end: u64,
        t_start: u64,
        t_end: u64,
        score: f64,
    }

    impl ChainItem for TestItem {
        fn q_start(&self) -> u64 { self.q_start }
        fn q_end(&self) -> u64 { self.q_end }
        fn t_start(&self) -> u64 { self.t_start }
        fn t_end(&self) -> u64 { self.t_end }
        fn score(&self) -> f64 { self.score }
    }

    #[test]
    fn test_kd_tree_build_and_search() {
        let items = vec![
            TestItem { q_start: 0, q_end: 10, t_start: 0, t_end: 10, score: 100.0 }, // 0
            TestItem { q_start: 20, q_end: 30, t_start: 20, t_end: 30, score: 100.0 }, // 1
            TestItem { q_start: 50, q_end: 60, t_start: 50, t_end: 60, score: 100.0 }, // 2
        ];

        let mut indices: Vec<usize> = (0..items.len()).collect();
        let mut tree = KdTree::build(&mut indices, &items);

        // Update scores in tree (simulating DP)
        // Assume item 0 has total score 100
        tree.update_scores(0, 100.0, &items);

        // Now search for predecessor for item 1
        // Cost func: simple distance penalty
        let cost_func = |cand_idx: usize, target_idx: usize| -> Option<f64> {
            if cand_idx >= target_idx { return None; } // strict order for test
            let cand = &items[cand_idx];
            let target = &items[target_idx];
            if cand.q_end > target.q_start || cand.t_end > target.t_start {
                return None;
            }
            let dist = (target.q_start - cand.q_end) + (target.t_start - cand.t_end);
            Some(100.0 + target.score - dist as f64) // 100.0 is prev total score
        };
        let lower_bound_func = |dq: u64, dt: u64| -> f64 {
            (dq + dt) as f64
        };

        let (best_score, best_pred) = tree.best_predecessor(
            1, // target is item 1
            100.0, // base score of item 1
            &items,
            &cost_func,
            &lower_bound_func
        );

        // Dist between 0 and 1:
        // 0 ends at 10, 10
        // 1 starts at 20, 20
        // dist = (20-10) + (20-10) = 20
        // score = 100 (prev) + 100 (curr) - 20 = 180
        assert_eq!(best_pred, Some(0));
        assert_eq!(best_score, 180.0);
    }
}
