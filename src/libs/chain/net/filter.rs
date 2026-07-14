//! Net filtering: prune net tree by score/size/type/synteny criteria.
//!
//! Mirrors UCSC `netFilter` semantics. A `FilterCriteria` holds optional
//! thresholds; `prune_gap` walks the tree and drops fills/gaps that fail.

use super::types::{Chrom, Fill, Gap};
use std::cell::{Ref, RefCell};
use std::collections::HashSet;
use std::rc::Rc;

/// Filtering criteria for net pruning (UCSC `netFilter` equivalent).
pub struct FilterCriteria {
    /// Minimum fill score.
    pub min_score: Option<f64>,
    /// Maximum fill score.
    pub max_score: Option<f64>,
    /// Minimum gap size to retain.
    pub min_gap: Option<u64>,
    /// Minimum aligned bases in a fill.
    pub min_ali: Option<u64>,
    /// Maximum aligned bases in a fill.
    pub max_ali: Option<u64>,
    /// Minimum target-side fill size.
    pub min_size_t: Option<u64>,
    /// Minimum query-side fill size.
    pub min_size_q: Option<u64>,
    /// Allowed target chromosome names.
    pub t_names: Option<HashSet<String>>,
    /// Excluded target chromosome names.
    pub not_t_names: Option<HashSet<String>>,
    /// Allowed query chromosome names.
    pub q_names: Option<HashSet<String>>,
    /// Excluded query chromosome names.
    pub not_q_names: Option<HashSet<String>>,
    /// Allowed synteny classes.
    pub types: Option<HashSet<String>>,

    // Synteny specific
    /// Keep only syntenic fills.
    pub do_syn: bool,
    /// Keep only non-syntenic fills.
    pub do_nonsyn: bool,
    /// Minimum score for a top-level fill to be considered syntenic.
    pub min_top_score: f64,
    /// Minimum score for a non-top fill to be considered syntenic.
    pub min_syn_score: f64,
    /// Minimum target size for a non-top fill to be considered syntenic.
    pub min_syn_size: f64,
    /// Minimum aligned bases for a non-top fill to be considered syntenic.
    pub min_syn_ali: u64,
    /// Maximum allowed qFar for syntenic fills.
    pub max_far: i64,

    /// Output only fill lines.
    pub fill_only: bool,
    /// Output only gap lines.
    pub gap_only: bool,
}

impl Default for FilterCriteria {
    fn default() -> Self {
        Self {
            min_score: None,
            max_score: None,
            min_gap: None,
            min_ali: None,
            max_ali: None,
            min_size_t: None,
            min_size_q: None,
            t_names: None,
            not_t_names: None,
            q_names: None,
            not_q_names: None,
            types: None,
            do_syn: false,
            do_nonsyn: false,
            min_top_score: 300000.0,
            min_syn_score: 200000.0,
            min_syn_size: 20000.0,
            min_syn_ali: 10000,
            max_far: 200000,
            fill_only: false,
            gap_only: false,
        }
    }
}

/// Check if a chromosome passes the target-name filter.
pub fn filter_chrom(chrom: &Chrom, c: &FilterCriteria) -> bool {
    if let Some(names) = &c.t_names {
        if !names.contains(&chrom.name) {
            return false;
        }
    }
    if let Some(names) = &c.not_t_names {
        if names.contains(&chrom.name) {
            return false;
        }
    }
    true
}

fn syn_filter(fill: &Fill, c: &FilterCriteria) -> bool {
    if fill.class.is_empty() {
        return false;
    }
    let t_size = fill.end - fill.start;

    if fill.score >= c.min_syn_score
        && (t_size as f64) >= c.min_syn_size
        && fill.ali >= c.min_syn_ali
    {
        return true;
    }
    if fill.class == "top" {
        return fill.score >= c.min_top_score;
    }
    if fill.class == "nonSyn" {
        return false;
    }
    if fill.q_far.unwrap_or(0) > c.max_far {
        return false;
    }
    true
}

fn filter_one(fill: &Fill, c: &FilterCriteria) -> bool {
    if let Some(names) = &c.q_names {
        if !names.contains(&fill.o_chrom) {
            return false;
        }
    }
    if let Some(names) = &c.not_q_names {
        if names.contains(&fill.o_chrom) {
            return false;
        }
    }
    if let Some(types) = &c.types {
        if !types.contains(&fill.class) {
            return false;
        }
    }

    if c.gap_only {
        return false;
    }

    if let Some(min_q) = c.min_size_q {
        let q_size = fill.o_end - fill.o_start;
        if q_size < min_q {
            return false;
        }
    }
    if let Some(min_t) = c.min_size_t {
        let t_size = fill.end - fill.start;
        if t_size < min_t {
            return false;
        }
    }

    if let Some(min_s) = c.min_score {
        if fill.score < min_s {
            return false;
        }
    }
    if let Some(max_s) = c.max_score {
        if fill.score > max_s {
            return false;
        }
    }

    if let Some(min_a) = c.min_ali {
        if fill.ali < min_a {
            return false;
        }
    }
    if let Some(max_a) = c.max_ali {
        if fill.ali > max_a {
            return false;
        }
    }

    if c.do_syn && !syn_filter(fill, c) {
        return false;
    }
    if c.do_nonsyn && syn_filter(fill, c) {
        return false;
    }

    true
}

/// Recursively prune a gap's fills (and their nested gaps) against the criteria.
pub fn prune_gap(gap: &Rc<RefCell<Gap>>, c: &FilterCriteria) {
    let mut gap_mut = gap.borrow_mut();

    let mut new_fills = Vec::new();

    for fill_rc in &gap_mut.fills {
        let keep = {
            let fill: Ref<Fill> = fill_rc.borrow();
            filter_one(&fill, c)
        };

        if keep {
            prune_fill(fill_rc, c);
            new_fills.push(fill_rc.clone());
        }
    }

    gap_mut.fills = new_fills;
}

fn prune_fill(fill: &Rc<RefCell<Fill>>, c: &FilterCriteria) {
    let mut fill_mut = fill.borrow_mut();

    let mut new_gaps = Vec::new();
    for gap_rc in &fill_mut.gaps {
        let keep = {
            let gap: Ref<Gap> = gap_rc.borrow();
            if c.fill_only {
                false
            } else if let Some(min_g) = c.min_gap {
                (gap.end - gap.start) >= min_g
            } else {
                true
            }
        };

        if keep {
            prune_gap(gap_rc, c);
            new_gaps.push(gap_rc.clone());
        }
    }
    fill_mut.gaps = new_gaps;
}
