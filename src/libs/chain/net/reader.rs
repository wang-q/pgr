//! UCSC Net text format reader.

use super::types::{Chrom, Fill, Gap, NetNode};
use std::cell::RefCell;
use std::io::{self, BufRead};
use std::rc::Rc;

pub fn read_nets<R: BufRead>(mut reader: R) -> io::Result<Vec<Chrom>> {
    let mut chroms = Vec::new();
    let mut current_chrom: Option<Chrom> = None;
    let mut stack: Vec<(usize, NetNode)> = Vec::new();
    let mut pending_comments = Vec::new();

    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        if line.trim().is_empty() {
            line.clear();
            continue;
        }

        if line.starts_with('#') {
            pending_comments.push(line.trim_end().to_string());
            line.clear();
            continue;
        }

        let mut indent = 0;
        for c in line.chars() {
            if c == ' ' {
                indent += 1;
            } else {
                break;
            }
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.is_empty() {
            line.clear();
            continue;
        }

        match parts[0] {
            "net" => {
                if let Some(c) = current_chrom {
                    chroms.push(c);
                }
                let name = parts[1];
                let size = parts[2].parse::<u64>().unwrap();
                let mut chrom = Chrom::new(name, size);
                if !pending_comments.is_empty() {
                    chrom.comments = std::mem::take(&mut pending_comments);
                }
                stack.clear();
                stack.push((0, NetNode::Gap(chrom.root.clone())));
                current_chrom = Some(chrom);
            }
            "fill" => {
                // fill tStart tLength qName qStrand qStart qLength id chainId score ali [type class]
                let start = parts[1].parse::<u64>().unwrap();
                let len = parts[2].parse::<u64>().unwrap();
                let q_name = parts[3].to_string();
                let q_strand = parts[4].chars().next().unwrap();
                let q_start = parts[5].parse::<u64>().unwrap();
                let q_len = parts[6].parse::<u64>().unwrap();
                // parts[7] is "id"
                let chain_id = parts[8].parse::<u64>().unwrap();
                // parts[9] is "score"
                let score = parts[10].parse::<f64>().unwrap();
                // parts[11] is "ali"
                let ali = parts[12].parse::<u64>().unwrap();

                let mut class = String::new();
                let mut q_dup = None;
                let mut q_over = None;
                let mut q_far = None;
                let mut t_n = None;
                let mut q_n = None;
                let mut t_r = None;
                let mut q_r = None;
                let mut t_trf = None;
                let mut q_trf = None;

                let mut i = 13;
                while i < parts.len() {
                    match parts[i] {
                        "type" => {
                            if i + 1 < parts.len() {
                                class = parts[i + 1].to_string();
                                i += 2;
                            } else {
                                i += 1;
                            }
                        }
                        "qDup" => {
                            let (v, ni) = parse_opt_u64(&parts, i);
                            q_dup = v;
                            i = ni;
                        }
                        "qOver" => {
                            let (v, ni) = parse_opt_u64(&parts, i);
                            q_over = v;
                            i = ni;
                        }
                        "qFar" => {
                            let (v, ni) = parse_opt_i64(&parts, i);
                            q_far = v;
                            i = ni;
                        }
                        "tN" => {
                            let (v, ni) = parse_opt_u64(&parts, i);
                            t_n = v;
                            i = ni;
                        }
                        "qN" => {
                            let (v, ni) = parse_opt_u64(&parts, i);
                            q_n = v;
                            i = ni;
                        }
                        "tR" => {
                            let (v, ni) = parse_opt_u64(&parts, i);
                            t_r = v;
                            i = ni;
                        }
                        "qR" => {
                            let (v, ni) = parse_opt_u64(&parts, i);
                            q_r = v;
                            i = ni;
                        }
                        "tTrf" => {
                            let (v, ni) = parse_opt_u64(&parts, i);
                            t_trf = v;
                            i = ni;
                        }
                        "qTrf" => {
                            let (v, ni) = parse_opt_u64(&parts, i);
                            q_trf = v;
                            i = ni;
                        }
                        _ => {
                            i += 1;
                        }
                    }
                }

                let fill = Rc::new(RefCell::new(Fill {
                    start,
                    end: start + len,
                    o_start: q_start,
                    o_end: q_start + q_len,
                    o_chrom: q_name,
                    o_strand: q_strand,
                    chain_id,
                    score,
                    ali,
                    class,
                    q_dup,
                    q_over,
                    q_far,
                    chain: None,
                    gaps: Vec::new(),
                    t_n,
                    q_n,
                    t_r,
                    q_r,
                    t_trf,
                    q_trf,
                }));

                // Find parent gap
                while let Some((parent_indent, parent_node)) = stack.last() {
                    if indent > *parent_indent {
                        if let NetNode::Gap(gap) = parent_node {
                            gap.borrow_mut().fills.push(fill.clone());
                            stack.push((indent, NetNode::Fill(fill)));
                            break;
                        } else {
                            stack.pop();
                        }
                    } else {
                        stack.pop();
                    }
                }
            }
            "gap" => {
                // gap tStart tLength qName qStrand qStart qLength
                let start = parts[1].parse::<u64>().unwrap();
                let len = parts[2].parse::<u64>().unwrap();
                let _q_name = parts[3].to_string();
                let _q_strand = parts[4].chars().next().unwrap();
                let q_start = parts[5].parse::<u64>().unwrap();
                let q_len = parts[6].parse::<u64>().unwrap();

                let mut t_n = None;
                let mut q_n = None;
                let mut t_r = None;
                let mut q_r = None;
                let mut t_trf = None;
                let mut q_trf = None;

                let mut i = 7;
                while i < parts.len() {
                    match parts[i] {
                        "tN" => {
                            let (v, ni) = parse_opt_u64(&parts, i);
                            t_n = v;
                            i = ni;
                        }
                        "qN" => {
                            let (v, ni) = parse_opt_u64(&parts, i);
                            q_n = v;
                            i = ni;
                        }
                        "tR" => {
                            let (v, ni) = parse_opt_u64(&parts, i);
                            t_r = v;
                            i = ni;
                        }
                        "qR" => {
                            let (v, ni) = parse_opt_u64(&parts, i);
                            q_r = v;
                            i = ni;
                        }
                        "tTrf" => {
                            let (v, ni) = parse_opt_u64(&parts, i);
                            t_trf = v;
                            i = ni;
                        }
                        "qTrf" => {
                            let (v, ni) = parse_opt_u64(&parts, i);
                            q_trf = v;
                            i = ni;
                        }
                        _ => {
                            i += 1;
                        }
                    }
                }

                let gap = Rc::new(RefCell::new(Gap {
                    start,
                    end: start + len,
                    o_start: q_start,
                    o_end: q_start + q_len,
                    fills: Vec::new(),
                    t_n,
                    q_n,
                    t_r,
                    q_r,
                    t_trf,
                    q_trf,
                }));

                // Find parent fill
                while let Some((parent_indent, parent_node)) = stack.last() {
                    if indent > *parent_indent {
                        if let NetNode::Fill(fill) = parent_node {
                            fill.borrow_mut().gaps.push(gap.clone());
                            stack.push((indent, NetNode::Gap(gap)));
                            break;
                        } else {
                            stack.pop();
                        }
                    } else {
                        stack.pop();
                    }
                }
            }
            _ => {}
        }
        line.clear();
    }
    if let Some(c) = current_chrom {
        chroms.push(c);
    }
    Ok(chroms)
}

// Parse an optional `name <value>` pair at position `i` in `parts` as u64.
fn parse_opt_u64(parts: &[&str], i: usize) -> (Option<u64>, usize) {
    if i + 1 < parts.len() {
        (Some(parts[i + 1].parse::<u64>().unwrap_or(0)), i + 2)
    } else {
        (None, i + 1)
    }
}

// Parse an optional `name <value>` pair at position `i` in `parts` as i64.
fn parse_opt_i64(parts: &[&str], i: usize) -> (Option<i64>, usize) {
    if i + 1 < parts.len() {
        (Some(parts[i + 1].parse::<i64>().unwrap_or(0)), i + 2)
    } else {
        (None, i + 1)
    }
}
