//! System resource helpers (memory, CPUs) backed by `sysinfo`.

use anyhow::{bail, Result};

/// Default in-memory sort budget (`--mem` default, 2 GiB).
pub const DEFAULT_MEM: u64 = 2 << 30;
/// Fraction of physical memory used as a hard safety ceiling.
pub const PHYS_MEM_FRACTION: f64 = 0.5;

/// Total physical memory in bytes, or `None` when unavailable.
pub fn physical_memory() -> Option<u64> {
    let sys = sysinfo::System::new_all();
    let bytes = sys.total_memory();
    (bytes > 0).then_some(bytes)
}

/// Number of logical CPUs (for future automatic parallelism).
pub fn logical_cpus() -> usize {
    let sys = sysinfo::System::new_all();
    sys.cpus().len()
}

/// Parses a BBTools-style KMG size ("2g", "512m", "1024", "2gb", "2gib")
/// into bytes; bare numbers are bytes. `k/m/g/t` are binary multiples.
pub fn parse_mem_size(s: &str) -> Result<u64> {
    let t = s.trim();
    if t.is_empty() {
        bail!("empty memory size");
    }
    let t = t.to_ascii_lowercase();
    // Strip an optional trailing "ib"/"b" first ("2gib", "2gb").
    let t = t
        .strip_suffix("ib")
        .or_else(|| t.strip_suffix('b'))
        .unwrap_or(&t);
    let (num, mult) = match t.as_bytes().last() {
        Some(b'k') => (&t[..t.len() - 1], 1u64 << 10),
        Some(b'm') => (&t[..t.len() - 1], 1u64 << 20),
        Some(b'g') => (&t[..t.len() - 1], 1u64 << 30),
        Some(b't') => (&t[..t.len() - 1], 1u64 << 40),
        _ => (t, 1u64),
    };
    if num.is_empty() {
        bail!("invalid memory size: {}", s);
    }
    let value: u64 = num
        .parse()
        .map_err(|_| anyhow::anyhow!("invalid memory size: {}", s))?;
    Ok(value.saturating_mul(mult))
}

/// Hard memory ceiling: min of the user cap and half the physical memory.
/// Falls back to the user cap when physical memory is unknown.
pub fn mem_cap(mem: Option<u64>) -> u64 {
    mem_cap_with(mem, physical_memory().unwrap_or(0))
}

/// Ceiling given an explicit physical-memory value (testable).
fn mem_cap_with(mem: Option<u64>, physical: u64) -> u64 {
    let cap = mem.unwrap_or(DEFAULT_MEM);
    if physical > 0 {
        cap.min((physical as f64 * PHYS_MEM_FRACTION) as u64)
    } else {
        cap
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_mem_size_accepts_kmg_and_plain() {
        assert_eq!(parse_mem_size("1024").unwrap(), 1024);
        assert_eq!(parse_mem_size("2g").unwrap(), 2 << 30);
        assert_eq!(parse_mem_size("512m").unwrap(), 512 << 20);
        assert_eq!(parse_mem_size("1k").unwrap(), 1024);
        assert_eq!(parse_mem_size("2gb").unwrap(), 2 << 30);
        assert_eq!(parse_mem_size("2gib").unwrap(), 2 << 30);
        assert!(parse_mem_size("").is_err());
        assert!(parse_mem_size("abc").is_err());
    }

    #[test]
    fn mem_cap_takes_the_min_of_cap_and_physical() {
        // The physical term must bind on small machines.
        let cap = mem_cap_with(Some(2 << 30), 1 << 30);
        assert_eq!(cap, 1 << 29); // 1 GiB machine -> half.
                                  // On big machines the user cap dominates.
        let cap = mem_cap_with(Some(2 << 30), 128 << 30);
        assert_eq!(cap, 2 << 30);
        // Unknown physical -> user cap only.
        let cap = mem_cap_with(Some(4 << 30), 0);
        assert_eq!(cap, 4 << 30);
    }
}
