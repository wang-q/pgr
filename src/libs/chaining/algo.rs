// use std::cmp::Ordering;

pub trait ChainItem {
    fn q_start(&self) -> u64;
    fn q_end(&self) -> u64;
    fn t_start(&self) -> u64;
    fn t_end(&self) -> u64;
    fn score(&self) -> f64;
}

pub struct KdLeaf<T> {
    pub item: T,
    pub original_idx: usize,
    pub best_pred_idx: Option<usize>,
    pub total_score: f64,
    pub hit: bool,
}

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

pub struct KdTree {
    root: Option<Box<KdNode>>,
}

impl KdTree {
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

    /// Find best predecessor.
    /// `cost_func` returns the connection cost (gap cost + overlap penalty).
    /// If connection is invalid (e.g. out of order), `cost_func` should return None.
    /// Returns (best_score, best_pred_idx)
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
