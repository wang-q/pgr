use std::fs::File;
use std::io::{BufRead, BufReader};
use anyhow::Result;

pub struct ScoreMatrix {
    matrix: Vec<i32>,
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl Default for ScoreMatrix {
    fn default() -> Self {
        // Default HoxD70 or similar
        let mut m = vec![0; 256 * 256];
        let bases = b"ACGT";
        // Simple identity matrix for default
        for &b1 in bases {
            for &b2 in bases {
                let idx = (b1 as usize) * 256 + (b2 as usize);
                if b1 == b2 {
                    m[idx] = 100; // Match
                } else {
                    m[idx] = -100; // Mismatch
                }
                // Handle lowercase
                let l1 = b1.to_ascii_lowercase();
                let l2 = b2.to_ascii_lowercase();
                
                m[(l1 as usize) * 256 + (l2 as usize)] = m[idx];
                m[(l1 as usize) * 256 + (b2 as usize)] = m[idx];
                m[(b1 as usize) * 256 + (l2 as usize)] = m[idx];
            }
        }
        // Handle N
        for i in 0..256 {
            m[('N' as usize) * 256 + i] = -100;
            m[i * 256 + ('N' as usize)] = -100;
            m[('n' as usize) * 256 + i] = -100;
            m[i * 256 + ('n' as usize)] = -100;
        }
        
        ScoreMatrix {
            matrix: m,
            gap_open: 400,
            gap_extend: 30,
        }
    }
}

impl ScoreMatrix {
    pub fn get_score(&self, c1: char, c2: char) -> i32 {
        self.matrix[(c1 as usize) * 256 + (c2 as usize)]
    }

    pub fn from_file(path: &str) -> Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut matrix = vec![0; 256 * 256];
        let mut gap_open = 400;
        let mut gap_extend = 30;
        
        // Default to A, C, G, T if no header found
        let mut chars = vec!['A', 'C', 'G', 'T'];
        let mut matrix_rows_read = 0;
        
        for line in reader.lines() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            
            // Check for gap scores at end: O = 400, E = 30
            // Also handle separate lines like O=400
            let has_o = line.contains("O=") || (line.contains("O") && line.contains("="));
            let has_e = line.contains("E=") || (line.contains("E") && line.contains("="));
            
            if has_o || has_e {
                let parts: Vec<&str> = line.split(|c| c == ',' || c == ' ' || c == '=').filter(|s| !s.is_empty()).collect();
                for i in 0..parts.len() {
                     if parts[i] == "O" && i+1 < parts.len() {
                         if let Ok(v) = parts[i+1].parse::<i32>() {
                             gap_open = v;
                         }
                     }
                     if parts[i] == "E" && i+1 < parts.len() {
                         if let Ok(v) = parts[i+1].parse::<i32>() {
                             gap_extend = v;
                         }
                     }
                }
                // Don't continue, might be mixed with other things? 
                // Usually these are parameter lines.
                continue;
            }
            
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.is_empty() { continue; }
            
            // Check if this is a header line (all single letters ACGT)
            if matrix_rows_read == 0 && parts.iter().all(|s| s.len() == 1 && "ACGTN".contains(s.chars().next().unwrap_or('?'))) {
                chars = parts.iter().map(|s| s.chars().next().unwrap()).collect();
                // eprintln!("Debug: Found header: {:?}", chars);
                continue;
            }
            
            // Read matrix row
            if matrix_rows_read < chars.len() {
                let row_char = chars[matrix_rows_read];
                // Check if line starts with the row char
                let val_start = if parts.len() > chars.len() { 1 } else { 0 };
                
                // Extra check: if val_start is 1, parts[0] must match row_char
                if val_start == 1 && parts[0].chars().next().unwrap() != row_char {
                    // Mismatch, maybe not a matrix line
                    continue;
                }

                let mut row_ok = false;
                for j in 0..chars.len() {
                    if j + val_start < parts.len() {
                        if let Ok(val) = parts[j + val_start].parse::<i32>() {
                            let col_char = chars[j];
                            let idx = (row_char as usize) * 256 + (col_char as usize);
                            matrix[idx] = val;
                            
                            // Fill lower case too
                            let r_lower = row_char.to_ascii_lowercase();
                            let c_lower = col_char.to_ascii_lowercase();
                            
                            matrix[(r_lower as usize) * 256 + (c_lower as usize)] = val;
                            matrix[(row_char as usize) * 256 + (c_lower as usize)] = val;
                            matrix[(r_lower as usize) * 256 + (col_char as usize)] = val;
                            row_ok = true;
                        }
                    }
                }
                if row_ok {
                    matrix_rows_read += 1;
                    // eprintln!("Debug: Read row {}, values: {:?}", row_char, parts);
                }
            }
        }
        
        // eprintln!("Debug: Matrix loaded. gap_open={}, gap_extend={}", gap_open, gap_extend);
        Ok(ScoreMatrix {
            matrix,
            gap_open,
            gap_extend,
        })
    }
}
