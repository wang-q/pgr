//! UCSC Net text format writer (filtered output with subchain scoring).
//!
//! `write_net` emits fills/gaps passing `min_score`/`min_fill` thresholds,
//! recomputing per-fill subchain score/size from the attached chain data.

use super::types::{Chrom, Fill, Gap};
use crate::libs::alignment::coords::reverse_range;
use crate::libs::chain::record::Chain;
use std::cell::RefCell;
use std::io::{self, Write};
use std::rc::Rc;

/// Length of the intersection of `[start1, end1)` and `[start2, end2)`.
pub fn range_intersection(start1: u64, end1: u64, start2: u64, end2: u64) -> u64 {
    let s = start1.max(start2);
    let e = end1.min(end2);
    e.saturating_sub(s)
}

pub fn write_net<W: Write>(
    chrom: &Chrom,
    writer: &mut W,
    is_q: bool,
    min_score: f64,
    min_fill: u64,
) -> io::Result<()> {
    if chrom.root.borrow().fills.is_empty() {
        return Ok(());
    }
    for comment in &chrom.comments {
        writeln!(writer, "{}", comment)?;
    }
    writeln!(writer, "net {} {}", chrom.name, chrom.size)?;

    for fill in &chrom.root.borrow().fills {
        write_fill(fill, writer, 1, is_q, min_score, min_fill)?;
    }
    Ok(())
}

fn write_fill<W: Write>(
    fill: &Rc<RefCell<Fill>>,
    writer: &mut W,
    depth: usize,
    is_q: bool,
    min_score: f64,
    min_fill: u64,
) -> io::Result<()> {
    let f = fill.borrow();

    // Calculate subscore/subsize if chain is available, otherwise use stored
    let (sub_size, sub_score) = if let Some(chain) = &f.chain {
        subchain_info(chain, f.start, f.end, is_q)
    } else {
        (f.ali, f.score)
    };

    if sub_score >= min_score && sub_size >= min_fill {
        write_indent(writer, depth)?;

        let mut line = format!(
            "fill {} {} {} {} {} {} id {} score {:.0} ali {}",
            f.start,
            f.end - f.start,
            f.o_chrom,
            f.o_strand,
            f.o_start,
            f.o_end - f.o_start,
            f.chain_id,
            sub_score,
            sub_size
        );

        // Optional fields: qOver, qFar, qDup, type
        // The order corresponds to UCSC's chainNet.c cnFillWrite implementation.
        push_opt_u64(&mut line, "qOver", f.q_over);
        push_opt_i64(&mut line, "qFar", f.q_far);
        push_opt_u64(&mut line, "qDup", f.q_dup);
        if !f.class.is_empty() {
            line.push_str(" type ");
            line.push_str(&f.class);
        }
        push_opt_u64(&mut line, "tN", f.t_n);
        push_opt_u64(&mut line, "qN", f.q_n);
        push_opt_u64(&mut line, "tR", f.t_r);
        push_opt_u64(&mut line, "qR", f.q_r);
        push_opt_u64(&mut line, "tTrf", f.t_trf);
        push_opt_u64(&mut line, "qTrf", f.q_trf);

        writeln!(writer, "{}", line)?;

        for gap in &f.gaps {
            write_gap(gap, fill, writer, depth + 1, is_q, min_score, min_fill)?;
        }
    }
    Ok(())
}

fn write_gap<W: Write>(
    gap: &Rc<RefCell<Gap>>,
    parent: &Rc<RefCell<Fill>>,
    writer: &mut W,
    depth: usize,
    is_q: bool,
    min_score: f64,
    min_fill: u64,
) -> io::Result<()> {
    let g = gap.borrow();
    let p = parent.borrow();
    let o_chrom = &p.o_chrom;
    let o_strand = p.o_strand;

    let mut line = format!(
        "gap {} {} {} {} {} {}",
        g.start,
        g.end - g.start,
        o_chrom,
        o_strand,
        g.o_start,
        g.o_end - g.o_start
    );
    push_opt_u64(&mut line, "tN", g.t_n);
    push_opt_u64(&mut line, "qN", g.q_n);
    push_opt_u64(&mut line, "tR", g.t_r);
    push_opt_u64(&mut line, "qR", g.q_r);
    push_opt_u64(&mut line, "tTrf", g.t_trf);
    push_opt_u64(&mut line, "qTrf", g.q_trf);

    write_indent(writer, depth)?;
    writeln!(writer, "{}", line)?;

    for fill in &g.fills {
        write_fill(fill, writer, depth + 1, is_q, min_score, min_fill)?;
    }
    Ok(())
}

fn write_indent<W: Write>(writer: &mut W, depth: usize) -> io::Result<()> {
    for _ in 0..depth {
        write!(writer, " ")?;
    }
    Ok(())
}

// Append `" name value"` to `line` if `val` is Some.
fn push_opt_u64(line: &mut String, name: &str, val: Option<u64>) {
    if let Some(v) = val {
        line.push(' ');
        line.push_str(name);
        line.push(' ');
        line.push_str(&v.to_string());
    }
}

// Append `" name value"` to `line` if `val` is Some.
fn push_opt_i64(line: &mut String, name: &str, val: Option<i64>) {
    if let Some(v) = val {
        line.push(' ');
        line.push_str(name);
        line.push(' ');
        line.push_str(&v.to_string());
    }
}

fn chain_base_count_sub_t(chain: &Chain, t_min: u64, t_max: u64) -> u64 {
    let mut total = 0;
    let mut t_curr = chain.header.t_start;
    for d in &chain.data {
        let t_start = t_curr;
        let t_end = t_curr + d.size;
        total += range_intersection(t_start, t_end, t_min, t_max);
        t_curr += d.size + d.dt;
    }
    total
}

fn chain_base_count_sub_q(chain: &Chain, q_min: u64, q_max: u64) -> u64 {
    let mut total = 0;
    let mut q_curr = chain.header.q_start;
    for d in &chain.data {
        let q_start = q_curr;
        let q_end = q_curr + d.size;
        total += range_intersection(q_start, q_end, q_min, q_max);
        q_curr += d.size + d.dq;
    }
    total
}

fn subchain_info(chain: &Chain, start: u64, end: u64, is_q: bool) -> (u64, f64) {
    let mut full_ali_size = 0;
    for d in &chain.data {
        full_ali_size += d.size;
    }

    if full_ali_size == 0 {
        return (0, 0.0);
    }

    let sub_size = if is_q {
        let mut s = start;
        let mut e = end;
        if chain.header.q_strand == '-' {
            reverse_range(&mut s, &mut e, chain.header.q_size);
        }
        if s <= chain.header.q_start && e >= chain.header.q_end {
            full_ali_size
        } else {
            chain_base_count_sub_q(chain, s, e)
        }
    } else if start <= chain.header.t_start && end >= chain.header.t_end {
        full_ali_size
    } else {
        chain_base_count_sub_t(chain, start, end)
    };

    let sub_score = chain.header.score * (sub_size as f64) / (full_ali_size as f64);
    (sub_size, sub_score)
}

/// Sort chroms by name and write each to `writer` via `finalize_net` + `write_net`.
pub fn write_sorted_net<W: Write>(
    net: &super::builder::ChainNet,
    writer: &mut W,
    is_q: bool,
    min_score: f64,
    min_fill: u64,
) -> anyhow::Result<()> {
    let mut chrom_names: Vec<_> = net.chroms.keys().cloned().collect();
    chrom_names.sort();
    for name in chrom_names {
        if let Some(chrom_cell) = net.chroms.get(&name) {
            let mut chrom = chrom_cell.borrow_mut();
            super::finalize::finalize_net(&mut chrom, is_q);
            write_net(&chrom, writer, is_q, min_score, min_fill)?;
        }
    }
    Ok(())
}

/// Write a net file with header comments and sorted net entries.
pub fn write_net_file(
    path: &str,
    net: &super::builder::ChainNet,
    is_q: bool,
    comments: &[String],
    min_score: f64,
    min_fill: u64,
) -> anyhow::Result<()> {
    use anyhow::Context;
    let mut writer = crate::libs::io::writer(path)
        .with_context(|| format!("Failed to open writer for {}", path))?;
    for comment in comments {
        write!(writer, "{}", comment)?;
    }
    write_sorted_net(net, &mut writer, is_q, min_score, min_fill)?;
    writer.flush()?;
    Ok(())
}
