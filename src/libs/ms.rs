use anyhow::Result;
use std::io::BufRead;
use std::process;
use std::time::{SystemTime, UNIX_EPOCH};

pub struct MsHeader {
    pub nsam: usize,
    pub howmany: usize,
    pub nsite: usize,
    pub npop: usize,
    pub sample_sizes: Option<Vec<usize>>,
}

pub struct MsSample {
    pub segsites: usize,
    pub positions: Vec<f64>,
    pub haplotypes: Vec<Vec<u8>>,
}

pub fn parse_header(line: &str) -> Result<MsHeader> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.is_empty() {
        anyhow::bail!("Empty header line");
    }
    let mut nsam = None;
    let mut howmany = None;
    for t in &tokens {
        if nsam.is_none() {
            if let Ok(v) = t.parse::<usize>() {
                nsam = Some(v);
                continue;
            }
        } else if howmany.is_none() {
            if let Ok(v) = t.parse::<usize>() {
                howmany = Some(v);
                break;
            }
        }
    }
    let nsam = nsam.ok_or_else(|| anyhow::anyhow!("Cannot parse nsam"))?;
    let howmany = howmany.ok_or_else(|| anyhow::anyhow!("Cannot parse howmany"))?;

    let mut nsite = 0usize;
    let mut npop = 1usize;
    let mut sample_sizes: Option<Vec<usize>> = None;
    let mut i = 0usize;
    while i < tokens.len() {
        match tokens[i] {
            "-r" => {
                if i + 2 < tokens.len() {
                    if let Ok(v) = tokens[i + 2].parse::<usize>() {
                        nsite = v;
                    }
                    i += 3;
                    continue;
                }
            }
            "-I" => {
                if i + 1 < tokens.len() {
                    if let Ok(v) = tokens[i + 1].parse::<usize>() {
                        npop = v;
                        let mut sizes = Vec::with_capacity(npop);
                        for k in 0..npop {
                            if i + 2 + k < tokens.len() {
                                sizes.push(tokens[i + 2 + k].parse::<usize>().unwrap_or(0));
                            }
                        }
                        sample_sizes = Some(sizes);
                    }
                    i += 2 + npop;
                    continue;
                }
            }
            _ => {}
        }
        i += 1;
    }

    Ok(MsHeader {
        nsam,
        howmany,
        nsite,
        npop,
        sample_sizes,
    })
}

pub fn read_next_sample<R: BufRead>(reader: &mut R, nsam: usize) -> Result<Option<MsSample>> {
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            return Ok(None);
        }
        if line.starts_with("//") {
            break;
        }
    }
    let mut segsites = 0usize;
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if line.starts_with("segsites:") {
            segsites = line
                .split_whitespace()
                .nth(1)
                .unwrap_or("0")
                .parse()
                .unwrap_or(0);
            break;
        }
    }
    let mut positions: Vec<f64> = Vec::new();
    if segsites > 0 {
        while positions.len() < segsites {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            let mut iter = line.split_whitespace();
            if line.starts_with("positions:") {
                iter.next();
            }
            for tok in iter {
                if positions.len() >= segsites {
                    break;
                }
                if let Ok(v) = tok.parse::<f64>() {
                    positions.push(v);
                }
            }
        }
    }
    let mut haplotypes: Vec<Vec<u8>> = Vec::with_capacity(nsam);
    for _ in 0..nsam {
        line.clear();
        reader.read_line(&mut line)?;
        haplotypes.push(line.trim().as_bytes().to_vec());
    }
    Ok(Some(MsSample {
        segsites,
        positions,
        haplotypes,
    }))
}

pub struct SimpleRng {
    state: u64,
}
impl SimpleRng {
    pub fn new(seed: u64) -> Self {
        let init = if seed == 0 { 0x9e3779b97f4a7c15 } else { seed };
        Self { state: init }
    }
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_mul(6364136223846793005).wrapping_add(1);
        self.state
    }
    pub fn next_f64(&mut self) -> f64 {
        let x = self.next_u64() >> 11;
        (x as f64) * (1.0 / ((1u64 << 53) as f64))
    }
}

pub fn system_seed() -> u64 {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs();
    secs ^ (process::id() as u64)
}

pub fn perturb_positions(positions: &mut [f64], rng: &mut SimpleRng) {
    for p in positions.iter_mut() {
        let r1 = rng.next_f64();
        let mut r2 = rng.next_f64();
        r2 /= 10000.0;
        if r1 > 0.5 && *p - r2 >= 0.0 {
            *p -= r2;
        } else if *p + r2 <= 1.0 {
            *p += r2;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    #[test]
    fn test_parse_header_basic() {
        let line = "ms 10 1 -t 5 -r 0 1000 -I 2 5 5";
        let hdr = parse_header(line).unwrap();
        assert_eq!(hdr.nsam, 10);
        assert_eq!(hdr.howmany, 1);
        assert_eq!(hdr.nsite, 1000);
        assert_eq!(hdr.npop, 2);
        assert_eq!(hdr.sample_sizes.as_ref().unwrap(), &vec![5, 5]);
    }

    #[test]
    fn test_read_next_sample_simple() {
        let nsam = 3;
        let input = "\
ms 3 1 -r 1.0 20
//
segsites: 2
positions: 0.1000 0.5000
01
10
11
//
segsites: 1
positions: 0.7500
1
0
1
";
        let mut reader = BufReader::new(input.as_bytes());
        let mut header = String::new();
        reader.read_line(&mut header).unwrap();
        let _hdr = parse_header(&header).unwrap();

        let s1 = read_next_sample(&mut reader, nsam).unwrap().unwrap();
        assert_eq!(s1.segsites, 2);
        assert_eq!(s1.positions, vec![0.1000, 0.5000]);
        assert_eq!(s1.haplotypes.len(), nsam);
        assert_eq!(std::str::from_utf8(&s1.haplotypes[0]).unwrap(), "01");

        let s2 = read_next_sample(&mut reader, nsam).unwrap().unwrap();
        assert_eq!(s2.segsites, 1);
        assert_eq!(s2.positions, vec![0.7500]);
        assert_eq!(s2.haplotypes.len(), nsam);
        assert_eq!(std::str::from_utf8(&s2.haplotypes[1]).unwrap(), "0");

        let s3 = read_next_sample(&mut reader, nsam).unwrap();
        assert!(s3.is_none());
    }

    #[test]
    fn test_perturb_positions_bounds() {
        let mut rng = SimpleRng::new(123);
        let mut pos = vec![0.0, 0.5, 1.0];
        perturb_positions(&mut pos, &mut rng);
        assert!(pos[0] >= 0.0 && pos[0] <= 1.0);
        assert!(pos[1] >= 0.0 && pos[1] <= 1.0);
        assert!(pos[2] >= 0.0 && pos[2] <= 1.0);
    }

    #[test]
    fn test_read_next_sample_multiline_positions() {
        let nsam = 2;
        let input = "\
ms 2 1 -r 1.0 10
//
segsites: 3
positions: 0.1000 0.5000
0.9000
01
10
";
        let mut reader = BufReader::new(input.as_bytes());
        let mut header = String::new();
        reader.read_line(&mut header).unwrap();
        let _hdr = parse_header(&header).unwrap();
        let s = read_next_sample(&mut reader, nsam).unwrap().unwrap();
        assert_eq!(s.segsites, 3);
        assert_eq!(s.positions, vec![0.1000, 0.5000, 0.9000]);
        assert_eq!(s.haplotypes.len(), nsam);
    }
}
