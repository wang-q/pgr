//! LZ-diff V2 encoder/decoder, ported from AGC v3.2.3's CLZDiff_V2.
//!
//! LZ77 variant: builds a hash table over a 2-bit encoded reference, encodes
//! matches as (position_diff, length), literals as single bytes, and N-runs
//! as a special marker sequence. V2 adds: equal-sequences optimization,
//! match-to-end optimization, '!' back-reference literal, and get_code_skip1
//! fast scanning.

use anyhow::{anyhow, Result};

// Constants mirroring C++ CLZDiffBase (lz_diff.h)
const EMPTY_KEY32: u32 = !0u32;
const EMPTY_KEY16: u16 = !0u16;
const MAX_LOAD_FACTOR: f64 = 0.7;
const MAX_NO_TRIES: u32 = 64;
const INVALID_SYMBOL: u8 = 31;
const N_CODE: u8 = 4;
const N_RUN_STARTER_CODE: u8 = 30;
const MIN_NRUN_LEN: u32 = 4;
const HASHING_STEP: u32 = 4;

/// Hash table: 16-bit or 32-bit positions (auto-selected by reference length).
enum HashTable {
    Short(Vec<u16>),
    Long(Vec<u32>),
}

/// LZ-diff V2 encoder/decoder.
pub struct LzDiff {
    reference: Vec<u8>,
    ht: Option<HashTable>,
    min_match_len: u32,
    key_len: u32,
    key_mask: u64,
    ht_size: u64,
    ht_mask: u64,
    short_ht_ver: bool,
    index_ready: bool,
}

/// Encode an ASCII DNA byte to 2-bit value (A=0, C=1, G=2, T=3, N=4, other=31).
fn encode_base(c: u8) -> u8 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' | b'U' | b'u' => 3,
        _ if c.is_ascii_alphabetic() => N_CODE,
        _ => INVALID_SYMBOL,
    }
}

/// Forward match length between two slices, up to max_len.
fn compare_fwd(p: &[u8], q: &[u8], max_len: u32) -> u32 {
    let n = (max_len as usize).min(p.len()).min(q.len());
    let mut len = 0usize;
    while len < n && p[len] == q[len] {
        len += 1;
    }
    len as u32
}

fn encode_literal(c: u8, encoded: &mut Vec<u8>) {
    encoded.push(b'A' + c);
}

fn encode_nrun(len: u32, encoded: &mut Vec<u8>) {
    encoded.push(N_RUN_STARTER_CODE);
    append_int(encoded, (len - MIN_NRUN_LEN) as i64);
    encoded.push(N_CODE);
}

fn is_literal(c: u8) -> bool {
    (b'A'..=b'A' + 20).contains(&c) || c == b'!'
}

fn is_nrun(c: u8) -> bool {
    c == N_RUN_STARTER_CODE
}

/// Append ASCII decimal representation of x to `out`.
fn append_int(out: &mut Vec<u8>, x: i64) {
    if x == 0 {
        out.push(b'0');
        return;
    }
    let mut x = x;
    if x < 0 {
        out.push(b'-');
        x = -x;
    }
    let mut tmp = [0u8; 20];
    let mut p = tmp.len();
    while x > 0 {
        p -= 1;
        tmp[p] = b'0' + (x % 10) as u8;
        x /= 10;
    }
    out.extend_from_slice(&tmp[p..]);
}

/// Read ASCII decimal from `data` starting at `pos`, advancing `pos`.
fn read_int(data: &[u8], pos: &mut usize) -> Result<i64> {
    let mut is_neg = false;
    if *pos < data.len() && data[*pos] == b'-' {
        is_neg = true;
        *pos += 1;
    }
    let mut x: i64 = 0;
    while *pos < data.len() && data[*pos].is_ascii_digit() {
        x = x * 10 + (data[*pos] - b'0') as i64;
        *pos += 1;
    }
    Ok(if is_neg { -x } else { x })
}

/// Hash a u64 key using murmurhash3 (pgr's existing pattern).
fn hash_key(x: u64) -> u64 {
    murmurhash3::murmurhash3_x64_128(&x.to_le_bytes(), 42).0
}

impl LzDiff {
    /// Create a new LzDiff with the given minimum match length (default 18).
    pub fn new(min_match_len: u32) -> Self {
        let key_len = min_match_len.saturating_sub(HASHING_STEP) + 1;
        let key_mask = if 2 * key_len >= 64 {
            !0u64
        } else {
            (!0u64) >> (64 - 2 * key_len)
        };
        Self {
            reference: Vec::new(),
            ht: None,
            min_match_len,
            key_len,
            key_mask,
            ht_size: 0,
            ht_mask: 0,
            short_ht_ver: false,
            index_ready: false,
        }
    }

    /// Store the reference (2-bit encoded) and decide hash table version.
    pub fn prepare(&mut self, reference: &[u8]) {
        self.short_ht_ver = (reference.len() as u64) / (HASHING_STEP as u64) < 65535;
        self.reference = reference.iter().map(|&c| encode_base(c)).collect();
        self.reference
            .resize(self.reference.len() + self.key_len as usize, INVALID_SYMBOL);
        self.index_ready = false;
        self.ht = None;
    }

    /// Build the hash table over the stored reference. Required only for
    /// encode; decode does not need the hash table.
    pub fn prepare_index(&mut self) {
        if self.index_ready {
            return;
        }

        // Count valid positions in sparse mode (every HASHING_STEP).
        let mut no_prev_valid: u32 = 0;
        let mut cnt_mod: u32 = 0;
        let key_len_mod = self.key_len % HASHING_STEP;
        let mut ht_size: u64 = 0;
        for &c in &self.reference {
            if c < 4 {
                no_prev_valid += 1;
            } else {
                no_prev_valid = 0;
            }
            cnt_mod += 1;
            if cnt_mod == HASHING_STEP {
                cnt_mod = 0;
            }
            if cnt_mod == key_len_mod && no_prev_valid >= self.key_len {
                ht_size += 1;
            }
        }

        // Round up to power of 2, apply load factor.
        ht_size = ((ht_size as f64) / MAX_LOAD_FACTOR) as u64;
        ht_size = ht_size.next_power_of_two();
        if ht_size < 8 {
            ht_size = 8;
        }
        self.ht_size = ht_size;
        self.ht_mask = ht_size - 1;

        if self.short_ht_ver {
            let mut ht16 = vec![EMPTY_KEY16; ht_size as usize];
            self.make_index16(&mut ht16);
            self.ht = Some(HashTable::Short(ht16));
        } else {
            let mut ht32 = vec![EMPTY_KEY32; ht_size as usize];
            self.make_index32(&mut ht32);
            self.ht = Some(HashTable::Long(ht32));
        }
        self.index_ready = true;
    }

    fn make_index16(&self, ht: &mut [u16]) {
        let ref_size = self.reference.len() as u32;
        let mut i: u32 = 0;
        while i + self.key_len < ref_size {
            let x = self.get_code(&self.reference[i as usize..]);
            if x != !0u64 {
                let pos = hash_key(x) & self.ht_mask;
                for j in 0..MAX_NO_TRIES {
                    let idx = ((pos + j as u64) & self.ht_mask) as usize;
                    if ht[idx] == EMPTY_KEY16 {
                        ht[idx] = (i / HASHING_STEP) as u16;
                        break;
                    }
                }
            }
            i += HASHING_STEP;
        }
    }

    fn make_index32(&self, ht: &mut [u32]) {
        let ref_size = self.reference.len() as u32;
        let mut i: u32 = 0;
        while i + self.key_len < ref_size {
            let x = self.get_code(&self.reference[i as usize..]);
            if x != !0u64 {
                let pos = hash_key(x) & self.ht_mask;
                for j in 0..MAX_NO_TRIES {
                    let idx = ((pos + j as u64) & self.ht_mask) as usize;
                    if ht[idx] == EMPTY_KEY32 {
                        ht[idx] = i / HASHING_STEP;
                        break;
                    }
                }
            }
            i += HASHING_STEP;
        }
    }

    /// Extract key_len 2-bit bases from `s` as a u64. Returns !0u64 if any
    /// base is non-ACGT (code > 3).
    fn get_code(&self, s: &[u8]) -> u64 {
        let mut x: u64 = 0;
        for &b in s.iter().take(self.key_len as usize) {
            if b > 3 {
                return !0u64;
            }
            x = (x << 2) | (b as u64);
        }
        x
    }

    /// Fast sliding-window: shift previous code left by 2, mask, and add the
    /// new base at s[key_len - 1].
    fn get_code_skip1(&self, code: u64, s: &[u8]) -> u64 {
        let new_base = s[self.key_len as usize - 1];
        if new_base > 3 {
            return !0u64;
        }
        ((code << 2) & self.key_mask) | (new_base as u64)
    }

    /// Detect a run of N_CODE starting at `s[0]`. Returns 0 if fewer than 3.
    fn get_nrun_len(s: &[u8], max_len: u32) -> u32 {
        if s.len() < 3 || s[0] != N_CODE || s[1] != N_CODE || s[2] != N_CODE {
            return 0;
        }
        let mut len = 3usize;
        let limit = (max_len as usize).min(s.len());
        while len < limit && s[len] == N_CODE {
            len += 1;
        }
        len as u32
    }

    /// Search the hash table for the best match starting at text position `i`.
    #[allow(clippy::too_many_arguments)]
    fn find_best_match(
        &self,
        ht_pos: u64,
        text: &[u8],
        i: u32,
        max_len: u32,
        no_prev_literals: u32,
        match_pos: &mut u32,
        len_bck: &mut u32,
        len_fwd: &mut u32,
    ) -> bool {
        *len_fwd = 0;
        *len_bck = 0;
        let mut min_to_update = self.min_match_len;
        let ht = self.ht.as_ref().expect("index not prepared");
        let mut cur_ht_pos = ht_pos;

        for _ in 0..MAX_NO_TRIES {
            let h_pos = match ht {
                HashTable::Short(ht16) => {
                    let val = ht16[cur_ht_pos as usize];
                    if val == EMPTY_KEY16 {
                        break;
                    }
                    (val as u32) * HASHING_STEP
                }
                HashTable::Long(ht32) => {
                    let val = ht32[cur_ht_pos as usize];
                    if val == EMPTY_KEY32 {
                        break;
                    }
                    val * HASHING_STEP
                }
            };

            let text_slice = &text[i as usize..];
            let ref_slice = &self.reference[h_pos as usize..];
            let f_len = compare_fwd(text_slice, ref_slice, max_len);

            if f_len >= self.key_len {
                // Backward match: compare text[i-b_len-1] with reference[h_pos-b_len-1]
                let mut b_len: u32 = 0;
                let b_limit = no_prev_literals.min(h_pos);
                while b_len < b_limit && b_len < i {
                    let t_idx = i as usize - b_len as usize - 1;
                    let r_idx = h_pos as usize - b_len as usize - 1;
                    if text[t_idx] != self.reference[r_idx] {
                        break;
                    }
                    b_len += 1;
                }

                if b_len + f_len > min_to_update {
                    *len_bck = b_len;
                    *len_fwd = f_len;
                    *match_pos = h_pos;
                    min_to_update = b_len + f_len;
                }
            }

            cur_ht_pos = (cur_ht_pos + 1) & self.ht_mask;
        }

        *len_bck + *len_fwd >= self.min_match_len
    }

    /// V2 match encoding: <diff_pos>[,<len-min_match_len>].
    fn encode_match(&self, ref_pos: u32, len: u32, pred_pos: u32, encoded: &mut Vec<u8>) {
        let dif_pos = ref_pos as i64 - pred_pos as i64;
        append_int(encoded, dif_pos);
        if len != !0u32 {
            encoded.push(b',');
            append_int(encoded, (len - self.min_match_len) as i64);
        }
        encoded.push(b'.');
    }

    /// Encode `text` against the prepared reference into `encoded`.
    pub fn encode(&mut self, text: &[u8], encoded: &mut Vec<u8>) {
        if !self.index_ready {
            self.prepare_index();
        }
        let text_size = text.len() as u32;
        encoded.clear();

        // Equal sequences optimization
        if text_size == self.reference.len() as u32 - self.key_len {
            let expected: Vec<u8> = text.iter().map(|&c| encode_base(c)).collect();
            if expected == self.reference[..text_size as usize] {
                return;
            }
        }

        encoded.reserve(text.len() / 64);

        let mut i: u32 = 0;
        let mut pred_pos: u32 = 0;
        let mut no_prev_literals: u32 = 0;
        let mut x_prev: u64 = !0u64;

        while i + self.key_len < text_size {
            let x = if x_prev != !0u64 && no_prev_literals > 0 {
                self.get_code_skip1(x_prev, &text[i as usize..])
            } else {
                self.get_code(&text[i as usize..])
            };
            x_prev = x;

            if x == !0u64 {
                let nrun_len = Self::get_nrun_len(&text[i as usize..], text_size - i);
                if nrun_len >= MIN_NRUN_LEN {
                    encode_nrun(nrun_len, encoded);
                    i += nrun_len;
                    no_prev_literals = 0;
                } else {
                    encode_literal(encode_base(text[i as usize]), encoded);
                    i += 1;
                    pred_pos += 1;
                    no_prev_literals += 1;
                }
                continue;
            }

            let ht_pos = hash_key(x) & self.ht_mask;
            let mut match_pos: u32 = 0;
            let mut len_bck: u32 = 0;
            let mut len_fwd: u32 = 0;
            let max_len = text_size - i;

            if !self.find_best_match(
                ht_pos,
                text,
                i,
                max_len,
                no_prev_literals,
                &mut match_pos,
                &mut len_bck,
                &mut len_fwd,
            ) {
                encode_literal(encode_base(text[i as usize]), encoded);
                i += 1;
                pred_pos += 1;
                no_prev_literals += 1;
                continue;
            }

            // Backward replacement: pop len_bck literals
            if len_bck > 0 {
                for _ in 0..len_bck {
                    encoded.pop();
                }
                match_pos -= len_bck;
                pred_pos -= len_bck;
                i -= len_bck;
            }

            // '!' back-ref: replace matching literals
            if match_pos == pred_pos {
                let e_size = encoded.len();
                let limit = e_size.min(match_pos as usize);
                for k in 1..=limit {
                    let idx = e_size - k;
                    let c = encoded[idx];
                    if !c.is_ascii_uppercase() {
                        break;
                    }
                    if c - b'A' == self.reference[match_pos as usize - k] {
                        encoded[idx] = b'!';
                    }
                }
            }

            // Match-to-end optimization
            let total_len = len_bck + len_fwd;
            let is_end = i + total_len == text_size
                && match_pos + total_len == self.reference.len() as u32 - self.key_len;
            let len = if is_end { !0u32 } else { total_len };
            self.encode_match(match_pos, len, pred_pos, encoded);

            pred_pos = match_pos + total_len;
            i += total_len;
            no_prev_literals = 0;
        }

        // Tail literals
        while i < text_size {
            encode_literal(encode_base(text[i as usize]), encoded);
            i += 1;
        }
    }

    /// Decode `encoded` using the stored reference, producing 2-bit coded
    /// output in `decoded`.
    pub fn decode(&self, encoded: &[u8], decoded: &mut Vec<u8>) -> Result<()> {
        decoded.clear();
        let mut pos = 0usize;
        let mut pred_pos: u32 = 0;
        let ref_len = self.reference.len() as u32;

        while pos < encoded.len() {
            let c = encoded[pos];
            if is_literal(c) {
                let decoded_c = if c == b'!' {
                    pos += 1;
                    if pred_pos as usize >= self.reference.len() {
                        return Err(anyhow!("decode: pred_pos out of range for '!' literal"));
                    }
                    self.reference[pred_pos as usize]
                } else {
                    pos += 1;
                    c - b'A'
                };
                decoded.push(decoded_c);
                pred_pos += 1;
            } else if is_nrun(c) {
                pos += 1;
                let raw_len = read_int(encoded, &mut pos)?;
                if pos >= encoded.len() || encoded[pos] != N_CODE {
                    return Err(anyhow!("decode: malformed N-run (missing N_CODE suffix)"));
                }
                pos += 1;
                let len = (raw_len + MIN_NRUN_LEN as i64) as usize;
                decoded.resize(decoded.len() + len, N_CODE);
            } else {
                // Match
                let raw_pos = read_int(encoded, &mut pos)?;
                let ref_pos = (raw_pos + pred_pos as i64) as u32;
                let len = if pos < encoded.len() && encoded[pos] == b',' {
                    pos += 1;
                    let raw_len = read_int(encoded, &mut pos)?;
                    if pos >= encoded.len() || encoded[pos] != b'.' {
                        return Err(anyhow!("decode: malformed match (missing '.')"));
                    }
                    pos += 1;
                    (raw_len + self.min_match_len as i64) as u32
                } else {
                    if pos >= encoded.len() || encoded[pos] != b'.' {
                        return Err(anyhow!("decode: malformed match (missing '.')"));
                    }
                    pos += 1;
                    !0u32
                };
                let actual_len = if len == !0u32 {
                    ref_len.checked_sub(ref_pos).ok_or_else(|| {
                        anyhow!("decode: ref_pos {} > ref_len {}", ref_pos, ref_len)
                    })?
                } else {
                    len
                };
                let end = (ref_pos + actual_len) as usize;
                if end > self.reference.len() {
                    return Err(anyhow!(
                        "decode: match [{}, {}) out of range (ref_len {})",
                        ref_pos,
                        end,
                        self.reference.len()
                    ));
                }
                decoded.extend_from_slice(&self.reference[ref_pos as usize..end]);
                pred_pos = ref_pos + actual_len;
            }
        }
        Ok(())
    }

    /// Estimate the encoded size without generating output. Stops early if
    /// the estimate exceeds `bound`.
    pub fn estimate(&mut self, text: &[u8], bound: u32) -> usize {
        if !self.index_ready {
            self.prepare_index();
        }
        let text_size = text.len() as u32;
        let mut est_cost: usize = 0;

        // Equal sequences optimization
        if text_size == self.reference.len() as u32 - self.key_len {
            let expected: Vec<u8> = text.iter().map(|&c| encode_base(c)).collect();
            if expected == self.reference[..text_size as usize] {
                return 0;
            }
        }

        let mut i: u32 = 0;
        let mut pred_pos: u32 = 0;
        let mut no_prev_literals: u32 = 0;
        let mut x_prev: u64 = !0u64;

        while i + self.key_len < text_size {
            if est_cost > bound as usize {
                return est_cost;
            }

            let x = if x_prev != !0u64 && no_prev_literals > 0 {
                self.get_code_skip1(x_prev, &text[i as usize..])
            } else {
                self.get_code(&text[i as usize..])
            };
            x_prev = x;

            if x == !0u64 {
                let nrun_len = Self::get_nrun_len(&text[i as usize..], text_size - i);
                if nrun_len >= MIN_NRUN_LEN {
                    est_cost += self.cost_nrun(nrun_len - MIN_NRUN_LEN) as usize;
                    i += nrun_len;
                    no_prev_literals = 0;
                } else {
                    est_cost += 1;
                    i += 1;
                    pred_pos += 1;
                    no_prev_literals += 1;
                }
                continue;
            }

            let ht_pos = hash_key(x) & self.ht_mask;
            let mut match_pos: u32 = 0;
            let mut len_bck: u32 = 0;
            let mut len_fwd: u32 = 0;
            let max_len = text_size - i;

            if !self.find_best_match(
                ht_pos,
                text,
                i,
                max_len,
                no_prev_literals,
                &mut match_pos,
                &mut len_bck,
                &mut len_fwd,
            ) {
                est_cost += 1;
                i += 1;
                pred_pos += 1;
                no_prev_literals += 1;
                continue;
            }

            let total_len = len_bck + len_fwd;
            let is_end = i + total_len == text_size
                && match_pos + total_len == self.reference.len() as u32 - self.key_len;
            let len = if is_end { !0u32 } else { total_len };
            est_cost += self.cost_match(match_pos, len, pred_pos) as usize;

            pred_pos = match_pos + total_len;
            i += total_len;
            no_prev_literals = 0;
        }

        est_cost += (text_size - i) as usize;
        est_cost
    }

    fn int_len(&self, x: i32) -> u32 {
        if x >= 0 {
            self.uint_len(x as u32)
        } else {
            1 + self.uint_len((-x) as u32)
        }
    }

    fn uint_len(&self, x: u32) -> u32 {
        if x < 10 {
            1
        } else if x < 100 {
            2
        } else if x < 1000 {
            3
        } else if x < 10000 {
            4
        } else if x < 100000 {
            5
        } else if x < 1000000 {
            6
        } else if x < 10000000 {
            7
        } else {
            8
        }
    }

    fn cost_nrun(&self, x: u32) -> u32 {
        2 + self.uint_len(x)
    }

    fn cost_match(&self, ref_pos: u32, len: u32, pred_pos: u32) -> u32 {
        let dif_pos = ref_pos as i32 - pred_pos as i32;
        let mut r = self.int_len(dif_pos);
        if len != !0u32 {
            r += 1 + self.uint_len(len - self.min_match_len);
        }
        r + 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_dna(len: usize, seed: u64) -> Vec<u8> {
        use rand::rngs::StdRng;
        use rand::Rng;
        use rand::SeedableRng;
        let mut rng = StdRng::seed_from_u64(seed);
        (0..len)
            .map(|_| {
                let r = rng.random_range(0u8..4);
                match r {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                }
            })
            .collect()
    }

    fn random_dna_with_n(len: usize, seed: u64, n_freq: f64) -> Vec<u8> {
        use rand::rngs::StdRng;
        use rand::Rng;
        use rand::SeedableRng;
        let mut rng = StdRng::seed_from_u64(seed);
        (0..len)
            .map(|_| {
                if rng.random_range(0.0..1.0) < n_freq {
                    b'N'
                } else {
                    match rng.random_range(0u8..4) {
                        0 => b'A',
                        1 => b'C',
                        2 => b'G',
                        _ => b'T',
                    }
                }
            })
            .collect()
    }

    /// Assert that encode -> decode produces the same 2-bit coded sequence.
    fn assert_roundtrip(reference: &[u8], text: &[u8], min_match_len: u32) {
        let mut lz = LzDiff::new(min_match_len);
        lz.prepare(reference);
        lz.prepare_index();
        let mut encoded = Vec::new();
        lz.encode(text, &mut encoded);
        let mut decoded = Vec::new();
        lz.decode(&encoded, &mut decoded).expect("decode failed");
        let expected: Vec<u8> = text.iter().map(|&c| encode_base(c)).collect();
        assert_eq!(decoded, expected, "roundtrip mismatch");
    }

    #[test]
    fn test_encode_decode_roundtrip_acgt() {
        let reference = random_dna(2000, 42);
        let text = random_dna(2000, 100);
        assert_roundtrip(&reference, &text, 18);
    }

    #[test]
    fn test_encode_decode_roundtrip_identical() {
        let reference = random_dna(2000, 42);
        // Same sequence -> empty delta
        let mut lz = LzDiff::new(18);
        lz.prepare(&reference);
        lz.prepare_index();
        let mut encoded = Vec::new();
        lz.encode(&reference, &mut encoded);
        assert!(
            encoded.is_empty(),
            "identical sequence should produce empty delta"
        );
    }

    #[test]
    fn test_encode_decode_roundtrip_with_n() {
        let reference = random_dna_with_n(2000, 42, 0.05);
        let text = random_dna_with_n(2000, 100, 0.05);
        assert_roundtrip(&reference, &text, 18);
    }

    #[test]
    fn test_encode_decode_roundtrip_lowercase() {
        let mut reference = random_dna(2000, 42);
        for b in reference.iter_mut() {
            *b = b.to_ascii_lowercase();
        }
        let text = random_dna(2000, 100);
        assert_roundtrip(&reference, &text, 18);
    }

    #[test]
    fn test_encode_decode_empty_text() {
        let reference = random_dna(100, 42);
        let mut lz = LzDiff::new(18);
        lz.prepare(&reference);
        lz.prepare_index();
        let mut encoded = Vec::new();
        lz.encode(b"", &mut encoded);
        let mut decoded = Vec::new();
        lz.decode(&encoded, &mut decoded).expect("decode failed");
        assert!(decoded.is_empty());
    }

    #[test]
    fn test_encode_decode_single_base() {
        let reference = b"ACGTACGTACGTACGTACGT";
        let text = b"A";
        assert_roundtrip(reference, text, 18);
    }

    #[test]
    fn test_encode_decode_pure_n() {
        let reference = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let text = b"NNNNNNNNNN";
        assert_roundtrip(reference, text, 18);
    }

    #[test]
    fn test_decode_without_index() {
        // decode should work without prepare_index
        let reference = random_dna(2000, 42);
        let text = random_dna(2000, 100);
        let mut lz = LzDiff::new(18);
        lz.prepare(&reference);
        lz.prepare_index();
        let mut encoded = Vec::new();
        lz.encode(&text, &mut encoded);

        // New LzDiff with same reference but no index
        let mut lz2 = LzDiff::new(18);
        lz2.prepare(&reference);
        // Do NOT call prepare_index
        let mut decoded = Vec::new();
        lz2.decode(&encoded, &mut decoded).expect("decode failed");
        let expected: Vec<u8> = text.iter().map(|&c| encode_base(c)).collect();
        assert_eq!(decoded, expected);
    }

    #[test]
    fn test_min_match_len_variants() {
        let reference = random_dna(2000, 42);
        let text = random_dna(2000, 100);
        for &mml in &[15u32, 18, 21] {
            assert_roundtrip(&reference, &text, mml);
        }
    }

    #[test]
    fn test_encode_decode_long_sequence() {
        // > 64KB to trigger 32-bit hash table
        let reference = random_dna(100_000, 42);
        let text = random_dna(100_000, 100);
        assert_roundtrip(&reference, &text, 18);
    }

    #[test]
    fn test_encode_decode_with_variants() {
        // Simulate SNP + indel variants
        let reference = random_dna(2000, 42);
        let mut text = reference.clone();
        // SNPs
        text[100] = b'G';
        text[500] = b'C';
        // Insertion
        text.insert(1000, b'A');
        // Deletion
        text.remove(1500);
        assert_roundtrip(&reference, &text, 18);
    }

    #[test]
    fn test_estimate_matches_encode_size() {
        let reference = random_dna(2000, 42);
        let text = random_dna(2000, 100);
        let mut lz = LzDiff::new(18);
        lz.prepare(&reference);
        lz.prepare_index();
        let mut encoded = Vec::new();
        lz.encode(&text, &mut encoded);
        let est = lz.estimate(&text, !0u32);
        assert_eq!(
            est,
            encoded.len(),
            "estimate should match actual encoded size"
        );
    }

    #[test]
    fn test_estimate_bound_early_return() {
        let reference = random_dna(2000, 42);
        let text = random_dna(2000, 100);
        let mut lz = LzDiff::new(18);
        lz.prepare(&reference);
        lz.prepare_index();
        let est = lz.estimate(&text, 0);
        // With bound=0, should return early with a positive cost (text != reference)
        assert!(est > 0);
    }
}
