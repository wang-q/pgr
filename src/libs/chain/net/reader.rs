//! UCSC Net text format reader.

use super::types::{Chrom, Fill, Gap, NetNode};
use anyhow::{anyhow, bail, Result};
use std::cell::RefCell;
use std::io::BufRead;
use std::rc::Rc;

pub fn read_nets<R: BufRead>(mut reader: R) -> Result<Vec<Chrom>> {
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
                if parts.len() < 3 {
                    bail!("net line needs at least 3 fields: {}", line.trim_end());
                }
                if let Some(c) = current_chrom {
                    chroms.push(c);
                }
                let name = parts[1];
                let size = parse_u64(&parts, 2, "net size")?;
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
                if parts.len() < 11 {
                    bail!("fill line needs at least 11 fields: {}", line.trim_end());
                }
                let start = parse_u64(&parts, 1, "fill tStart")?;
                let len = parse_u64(&parts, 2, "fill tLength")?;
                let q_name = parts[3].to_string();
                let q_strand = parts[4]
                    .chars()
                    .next()
                    .ok_or_else(|| anyhow!("empty fill qStrand field"))?;
                let q_start = parse_u64(&parts, 5, "fill qStart")?;
                let q_len = parse_u64(&parts, 6, "fill qLength")?;
                // parts[7] is "id"
                let chain_id = parse_u64(&parts, 8, "fill chainId")?;
                // parts[9] is "score"
                let score = parse_f64(&parts, 10, "fill score")?;
                // parts[11] is "ali"
                let ali = parse_u64(&parts, 12, "fill ali")?;

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
                            let (v, ni) = parse_opt_u64(&parts, i)?;
                            q_dup = v;
                            i = ni;
                        }
                        "qOver" => {
                            let (v, ni) = parse_opt_u64(&parts, i)?;
                            q_over = v;
                            i = ni;
                        }
                        "qFar" => {
                            let (v, ni) = parse_opt_i64(&parts, i)?;
                            q_far = v;
                            i = ni;
                        }
                        "tN" => {
                            let (v, ni) = parse_opt_u64(&parts, i)?;
                            t_n = v;
                            i = ni;
                        }
                        "qN" => {
                            let (v, ni) = parse_opt_u64(&parts, i)?;
                            q_n = v;
                            i = ni;
                        }
                        "tR" => {
                            let (v, ni) = parse_opt_u64(&parts, i)?;
                            t_r = v;
                            i = ni;
                        }
                        "qR" => {
                            let (v, ni) = parse_opt_u64(&parts, i)?;
                            q_r = v;
                            i = ni;
                        }
                        "tTrf" => {
                            let (v, ni) = parse_opt_u64(&parts, i)?;
                            t_trf = v;
                            i = ni;
                        }
                        "qTrf" => {
                            let (v, ni) = parse_opt_u64(&parts, i)?;
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
                let mut parent_found = false;
                while let Some((parent_indent, parent_node)) = stack.last() {
                    if indent > *parent_indent {
                        if let NetNode::Gap(gap) = parent_node {
                            gap.borrow_mut().fills.push(fill.clone());
                            stack.push((indent, NetNode::Fill(fill)));
                            parent_found = true;
                            break;
                        } else {
                            stack.pop();
                        }
                    } else {
                        stack.pop();
                    }
                }
                if !parent_found {
                    bail!("orphaned fill line: {}", line.trim_end());
                }
            }
            "gap" => {
                // gap tStart tLength qName qStrand qStart qLength
                if parts.len() < 7 {
                    bail!("gap line needs at least 7 fields: {}", line.trim_end());
                }
                let start = parse_u64(&parts, 1, "gap tStart")?;
                let len = parse_u64(&parts, 2, "gap tLength")?;
                let _q_name = parts[3].to_string();
                let _q_strand = parts[4]
                    .chars()
                    .next()
                    .ok_or_else(|| anyhow!("empty gap qStrand field"))?;
                let q_start = parse_u64(&parts, 5, "gap qStart")?;
                let q_len = parse_u64(&parts, 6, "gap qLength")?;

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
                            let (v, ni) = parse_opt_u64(&parts, i)?;
                            t_n = v;
                            i = ni;
                        }
                        "qN" => {
                            let (v, ni) = parse_opt_u64(&parts, i)?;
                            q_n = v;
                            i = ni;
                        }
                        "tR" => {
                            let (v, ni) = parse_opt_u64(&parts, i)?;
                            t_r = v;
                            i = ni;
                        }
                        "qR" => {
                            let (v, ni) = parse_opt_u64(&parts, i)?;
                            q_r = v;
                            i = ni;
                        }
                        "tTrf" => {
                            let (v, ni) = parse_opt_u64(&parts, i)?;
                            t_trf = v;
                            i = ni;
                        }
                        "qTrf" => {
                            let (v, ni) = parse_opt_u64(&parts, i)?;
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
                let mut parent_found = false;
                while let Some((parent_indent, parent_node)) = stack.last() {
                    if indent > *parent_indent {
                        if let NetNode::Fill(fill) = parent_node {
                            fill.borrow_mut().gaps.push(gap.clone());
                            stack.push((indent, NetNode::Gap(gap)));
                            parent_found = true;
                            break;
                        } else {
                            stack.pop();
                        }
                    } else {
                        stack.pop();
                    }
                }
                if !parent_found {
                    bail!("orphaned gap line: {}", line.trim_end());
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

// Parse a required u64 field at index `i`.
fn parse_u64(parts: &[&str], i: usize, field: &str) -> Result<u64> {
    let s = parts
        .get(i)
        .copied()
        .ok_or_else(|| anyhow!("missing {field} at index {i}"))?;
    s.parse::<u64>()
        .map_err(|_| anyhow!("invalid {field} value: {s}"))
}

// Parse a required f64 field at index `i`.
fn parse_f64(parts: &[&str], i: usize, field: &str) -> Result<f64> {
    let s = parts
        .get(i)
        .copied()
        .ok_or_else(|| anyhow!("missing {field} at index {i}"))?;
    s.parse::<f64>()
        .map_err(|_| anyhow!("invalid {field} value: {s}"))
}

// Parse an optional `name <value>` pair at position `i` in `parts` as u64.
fn parse_opt_u64(parts: &[&str], i: usize) -> Result<(Option<u64>, usize)> {
    if i + 1 < parts.len() {
        let val = parts[i + 1].parse::<u64>().map_err(|e| {
            anyhow!(
                "invalid u64 value at index {}: {}: {}",
                i + 1,
                parts[i + 1],
                e
            )
        })?;
        Ok((Some(val), i + 2))
    } else {
        Ok((None, i + 1))
    }
}

// Parse an optional `name <value>` pair at position `i` in `parts` as i64.
fn parse_opt_i64(parts: &[&str], i: usize) -> Result<(Option<i64>, usize)> {
    if i + 1 < parts.len() {
        let val = parts[i + 1].parse::<i64>().map_err(|e| {
            anyhow!(
                "invalid i64 value at index {}: {}: {}",
                i + 1,
                parts[i + 1],
                e
            )
        })?;
        Ok((Some(val), i + 2))
    } else {
        Ok((None, i + 1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_malformed_net_missing_fields() {
        let data = "net chr1\n";
        let r = read_nets(std::io::Cursor::new(data));
        assert!(r.is_err());
    }

    #[test]
    fn test_malformed_net_non_numeric_size() {
        let data = "net chr1 notanumber\n";
        let r = read_nets(std::io::Cursor::new(data));
        assert!(r.is_err());
    }

    #[test]
    fn test_malformed_fill_missing_fields() {
        let data = "net chr1 100\nfill 10 20\n";
        let r = read_nets(std::io::Cursor::new(data));
        assert!(r.is_err());
    }

    #[test]
    fn test_malformed_fill_non_numeric() {
        let data = "net chr1 100\nfill abc 20 chr2 + 0 20 id 1 score 100 ali 20\n";
        let r = read_nets(std::io::Cursor::new(data));
        assert!(r.is_err());
    }

    #[test]
    fn test_malformed_gap_missing_fields() {
        let data = "net chr1 100\ngap 10 20 chr2\n";
        let r = read_nets(std::io::Cursor::new(data));
        assert!(r.is_err());
    }

    #[test]
    fn test_empty_input() {
        let r = read_nets(std::io::Cursor::new(""));
        assert!(r.is_ok());
        assert!(r.unwrap().is_empty());
    }

    #[test]
    fn test_binary_input_does_not_panic() {
        let binary = b"\xff\xfe\x00\x01net chr1 100\n";
        let _ = read_nets(std::io::Cursor::new(binary.as_slice()));
    }

    #[test]
    fn test_orphaned_fill() {
        let data = "net chr1 100\nfill 0 10 chr2 + 0 10 id 1 score 100 ali 10\n";
        let r = read_nets(std::io::Cursor::new(data));
        match r {
            Err(e) => assert!(e.to_string().contains("orphaned fill line")),
            Ok(_) => panic!("expected orphaned fill error"),
        }
    }

    #[test]
    fn test_orphaned_gap() {
        let data = "net chr1 100\ngap 0 10 chr2 + 0 10\n";
        let r = read_nets(std::io::Cursor::new(data));
        match r {
            Err(e) => assert!(e.to_string().contains("orphaned gap line")),
            Ok(_) => panic!("expected orphaned gap error"),
        }
    }

    #[test]
    fn test_orphaned_nested_fill() {
        // A gap at indent 1 has no parent fill at indent 0.
        let data = "net chr1 100\n gap 0 10 chr2 + 0 10\n";
        let r = read_nets(std::io::Cursor::new(data));
        assert!(r.is_err());
    }
}
