use super::align::{Alignment, AlignmentEngine, AlignmentParams, AlignmentType};
use super::graph::PoaGraph;
use petgraph::graph::NodeIndex;
use std::collections::HashMap;
use wide::i32x8;

/// Number of i32 lanes per SIMD vector (AVX2 256-bit).
const LANES: usize = 8;

/// Negative infinity sentinel, identical to the scalar engine.
const NEG_INF: i32 = -1_000_000_000;

/// Threshold below which values are treated as negative infinity (after
/// `NEG_INF + penalty` drift from direct vector arithmetic).
const NEG_INF_HALF: i32 = NEG_INF / 2;

/// Per-alignment context shared by both SIMD paths.
struct AlignCtx {
    sorted_nodes: Vec<NodeIndex>,
    node_map: HashMap<NodeIndex, usize>,
    preds_of_rank: Vec<Vec<usize>>,
    is_start: Vec<bool>,
    base_of_rank: Vec<u8>,
    n_nodes: usize,
    n_seq: usize,
    n_vec: usize,
}

/// SIMD partial-order alignment engine.
///
/// Platform policy (same as `libs::hv`): a hand-written AVX2 path is the
/// primary x86-64 implementation, dispatched at runtime with
/// `is_x86_feature_detected!`; all other targets fall through to the portable
/// `wide` implementation. Both paths are bit-identical to the scalar engine.
pub struct SimdAlignmentEngine {
    params: AlignmentParams,
    align_type: AlignmentType,
}

/// SIMD implementation selector, mainly for benchmarking and diagnostics.
pub enum SimdPath {
    /// Runtime detection: AVX2 on capable x86-64, `wide` elsewhere.
    Auto,
    /// Force the portable `wide` path.
    Wide,
    /// Force the AVX2 path (falls back to `wide` on CPUs without AVX2).
    Avx2,
}

impl SimdAlignmentEngine {
    /// Creates a new SIMD alignment engine.
    pub fn new(params: AlignmentParams, align_type: AlignmentType) -> Self {
        Self { params, align_type }
    }

    /// Aligns `sequence` to `graph` using the given implementation path.
    pub fn align_with(&self, path: SimdPath, sequence: &[u8], graph: &PoaGraph) -> Alignment {
        let ctx = build_ctx(sequence, graph);
        if ctx.n_nodes == 0 {
            return Alignment {
                score: 0,
                path: (0..ctx.n_seq).map(|i| (Some(i), None)).collect(),
            };
        }
        match path {
            SimdPath::Wide => align_wide(&self.params, self.align_type, sequence, graph, &ctx),
            SimdPath::Auto | SimdPath::Avx2 => {
                #[cfg(target_arch = "x86_64")]
                if is_x86_feature_detected!("avx2") {
                    // SAFETY: gated on runtime AVX2 support.
                    return unsafe {
                        align_avx2(&self.params, self.align_type, sequence, graph, &ctx)
                    };
                }
                align_wide(&self.params, self.align_type, sequence, graph, &ctx)
            }
        }
    }
}

impl AlignmentEngine for SimdAlignmentEngine {
    fn align(&self, sequence: &[u8], graph: &PoaGraph) -> Alignment {
        self.align_with(SimdPath::Auto, sequence, graph)
    }
}

fn build_ctx(sequence: &[u8], graph: &PoaGraph) -> AlignCtx {
    let sorted_nodes = graph.topological_sort();
    let node_map: HashMap<NodeIndex, usize> = sorted_nodes
        .iter()
        .enumerate()
        .map(|(i, &n)| (n, i))
        .collect();
    let preds_of_rank: Vec<Vec<usize>> = sorted_nodes
        .iter()
        .map(|&n| {
            graph
                .graph
                .neighbors_directed(n, petgraph::Direction::Incoming)
                .map(|p| node_map[&p])
                .collect()
        })
        .collect();
    let is_start: Vec<bool> = preds_of_rank.iter().map(|p| p.is_empty()).collect();
    let base_of_rank: Vec<u8> = sorted_nodes.iter().map(|&n| graph.graph[n].base).collect();
    let n_seq = sequence.len();
    let n_nodes = sorted_nodes.len();
    AlignCtx {
        sorted_nodes,
        node_map,
        preds_of_rank,
        is_start,
        base_of_rank,
        n_nodes,
        n_seq,
        n_vec: (n_seq + 1).div_ceil(LANES),
    }
}

/// Match/mismatch profile per vector column (`j == 0` and padding lanes are 0).
fn build_profile(
    seq: &[u8],
    base: u8,
    m_score: i32,
    n_score: i32,
    n_seq: usize,
    n_vec: usize,
) -> Vec<[i32; LANES]> {
    (0..n_vec)
        .map(|v| {
            let mut arr = [0i32; LANES];
            for (k, a) in arr.iter_mut().enumerate() {
                let j = v * LANES + k;
                *a = if j == 0 || j > n_seq {
                    0
                } else if seq[j - 1] == base {
                    m_score
                } else {
                    n_score
                };
            }
            arr
        })
        .collect()
}

/// Virtual-root row: `root[0] = 0`, `root[j] = gap_open + (j-1)*gap_extend`.
fn build_root(gap_open: i32, gap_extend: i32, n_seq: usize, n_vec: usize) -> Vec<[i32; LANES]> {
    (0..n_vec)
        .map(|v| {
            let mut arr = [0i32; LANES];
            for (k, a) in arr.iter_mut().enumerate() {
                let j = v * LANES + k;
                *a = if j == 0 || j > n_seq {
                    0
                } else {
                    gap_open + (j as i32 - 1) * gap_extend
                };
            }
            arr
        })
        .collect()
}

/// Shifts a vector right by one lane, filling lane 0 with `prev_last`:
/// `[prev_last, v0, v1, ..., v6]` (the `j-1` diagonal for a column).
#[inline]
fn shift1(v: i32x8, prev_last: i32) -> i32x8 {
    let arr = v.to_array();
    i32x8::from([
        prev_last, arr[0], arr[1], arr[2], arr[3], arr[4], arr[5], arr[6],
    ])
}

/// Prefix scan with linear penalties: `y[k] = max_{t <= k}(a[t] + (k-t)*e)`.
///
/// Log-step scan (shift 1/2/4 lanes), equivalent to spoa's `prefix_max`.
fn prefix_scan(mut a: i32x8, e: i32) -> i32x8 {
    let neg_inf = i32x8::splat(NEG_INF);
    for n in [1usize, 2, 4] {
        let arr = a.to_array();
        let mut shifted = [0i32; LANES];
        let mut excl = [0i32; LANES];
        for k in 0..LANES {
            shifted[k] = if k >= n { arr[k - n] } else { arr[LANES - 1] };
            excl[k] = if k >= n { 0 } else { -1 };
        }
        let excl_mask = i32x8::from(excl);
        let contrib = excl_mask.select(neg_inf, i32x8::from(shifted)) + i32x8::splat(n as i32 * e);
        a = a.max(contrib);
    }
    a
}

/// Fills lane 0 of `v` with `boundary`.
#[inline]
fn set_lane0(v: i32x8, boundary: i32) -> i32x8 {
    let arr = v.to_array();
    let mut out = arr;
    out[0] = boundary;
    i32x8::from(out)
}

/// Clamps negative-infinity-domain lanes to `NEG_INF`.
#[inline]
fn clamp_neg_inf(v: i32x8) -> i32x8 {
    let mask = (v - i32x8::splat(NEG_INF_HALF)).is_negative();
    mask.select(i32x8::splat(NEG_INF), v)
}

/// Forward DP (wide path); returns M/E/F matrices and best-cell info.
#[allow(clippy::type_complexity)]
fn forward_wide(
    params: &AlignmentParams,
    align_type: AlignmentType,
    sequence: &[u8],
    ctx: &AlignCtx,
) -> (
    Vec<Vec<i32x8>>,
    Vec<Vec<i32x8>>,
    Vec<Vec<i32x8>>,
    i32,
    usize,
    usize,
    u8,
) {
    let m_score = params.match_score;
    let n_score = params.mismatch_score;
    let g_open = params.gap_open;
    let g_ext = params.gap_extend;
    let is_local = align_type == AlignmentType::Local;
    let is_semi = align_type == AlignmentType::SemiGlobal;

    let n_nodes = ctx.n_nodes;
    let n_vec = ctx.n_vec;
    let g_open_v = i32x8::splat(g_open);
    let g_ext_v = i32x8::splat(g_ext);
    let neg_inf_v = i32x8::splat(NEG_INF);
    let j1_zero = i32x8::from([
        NEG_INF, 0, NEG_INF, NEG_INF, NEG_INF, NEG_INF, NEG_INF, NEG_INF,
    ]);

    let root: Vec<i32x8> = build_root(g_open, g_ext, ctx.n_seq, n_vec)
        .into_iter()
        .map(i32x8::from)
        .collect();

    let mut m = vec![vec![neg_inf_v; n_vec]; n_nodes];
    let mut e = vec![vec![neg_inf_v; n_vec]; n_nodes];
    let mut f = vec![vec![neg_inf_v; n_vec]; n_nodes];
    let mut f0_of = vec![NEG_INF; n_nodes];

    let mut best_score = if is_local { 0 } else { NEG_INF };
    let mut best_i = 0usize;
    let mut best_j = ctx.n_seq;
    let mut best_state = 0u8;

    for i in 0..n_nodes {
        let preds = &ctx.preds_of_rank[i];
        let is_start = ctx.is_start[i];
        let profile: Vec<i32x8> = build_profile(
            sequence,
            ctx.base_of_rank[i],
            m_score,
            n_score,
            ctx.n_seq,
            n_vec,
        )
        .into_iter()
        .map(i32x8::from)
        .collect();

        // Boundary value for F at j == 0 (column 0, lane 0).
        let f0 = if is_local || is_semi {
            0
        } else if is_start {
            g_open
        } else {
            let mut max_prev = NEG_INF;
            for &u in preds {
                max_prev = max_prev.max(f0_of[u]);
            }
            if max_prev > NEG_INF {
                max_prev + g_ext
            } else {
                NEG_INF
            }
        };
        f0_of[i] = f0;

        let mut m_last = NEG_INF;
        let mut e_last = NEG_INF;
        let mut f_last = NEG_INF;

        for v in 0..n_vec {
            let prof = profile[v];

            // M: match/mismatch from predecessors (diagonal j-1) or virtual root.
            let mut m_vec = neg_inf_v;
            if is_start {
                let root_last = if v == 0 {
                    NEG_INF
                } else {
                    root[v - 1].to_array()[LANES - 1]
                };
                m_vec = shift1(root[v], root_last) + prof;
            } else {
                for &u in preds {
                    let pm = m[u][v].max(e[u][v]).max(f[u][v]);
                    let last = if v == 0 {
                        NEG_INF
                    } else {
                        m[u][v - 1].max(e[u][v - 1]).max(f[u][v - 1]).to_array()[LANES - 1]
                    };
                    m_vec = m_vec.max(shift1(pm, last) + prof);
                }
                if (is_local || is_semi) && v == 0 {
                    m_vec = m_vec.max(j1_zero + prof);
                }
            }
            if is_local {
                m_vec = m_vec.max(i32x8::splat(0));
            }
            m_vec = clamp_neg_inf(m_vec);
            if v == 0 {
                m_vec = set_lane0(m_vec, NEG_INF);
            }

            // F: graph deletion from predecessor rows (same column j).
            let mut f_vec = neg_inf_v;
            if !is_start {
                for &u in preds {
                    let fm = m[u][v] + g_open_v;
                    let ff = f[u][v] + g_ext_v;
                    let fe = e[u][v] + g_open_v;
                    f_vec = f_vec.max(fm.max(ff).max(fe));
                }
            }
            if is_local {
                f_vec = f_vec.is_negative().select(neg_inf_v, f_vec);
            }
            f_vec = clamp_neg_inf(f_vec);
            if v == 0 {
                f_vec = set_lane0(f_vec, f0);
            }

            // E: sequence insertion; affine gap chain via prefix scan.
            let base_e = shift1(m_vec, m_last).max(shift1(f_vec, f_last)) + g_open_v;
            let x = if v == 0 {
                set_lane0(base_e, NEG_INF)
            } else {
                // Column start (lane 0, j = v*8) also continues E[j-1] + e
                // from the previous column's tail.
                let x0 = base_e.to_array()[0].max(e_last + g_ext);
                set_lane0(base_e, x0)
            };
            let mut e_vec = prefix_scan(x, g_ext);
            if is_local {
                e_vec = e_vec.is_negative().select(neg_inf_v, e_vec);
            }
            e_vec = clamp_neg_inf(e_vec);

            m[i][v] = m_vec;
            e[i][v] = e_vec;
            f[i][v] = f_vec;

            // Best score scan (local: all cells, last occurrence wins).
            if is_local {
                let col = m_vec.max(e_vec).max(f_vec);
                let col_arr = col.to_array();
                // Only consider valid lanes: j in 1..=n_seq.
                let start_k = if v == 0 { 1 } else { 0 };
                let end_k = (ctx.n_seq + 1 - v * LANES).min(LANES);
                if end_k > start_k {
                    let cm = col_arr[start_k..end_k].iter().copied().max().unwrap();
                    if cm >= best_score {
                        best_score = cm;
                        for k in (start_k..end_k).rev() {
                            if col_arr[k] == cm {
                                best_i = i;
                                best_j = v * LANES + k;
                                let m_arr = m_vec.to_array();
                                let e_arr = e_vec.to_array();
                                best_state = if m_arr[k] == cm {
                                    0
                                } else if e_arr[k] == cm {
                                    1
                                } else {
                                    2
                                };
                                break;
                            }
                        }
                    }
                }
            }

            m_last = m_vec.to_array()[LANES - 1];
            e_last = e_vec.to_array()[LANES - 1];
            f_last = f_vec.to_array()[LANES - 1];
        }
    }

    if !is_local {
        // Semi-global/global: last column (j = n_seq) only, first max wins.
        let v_last = ctx.n_seq / LANES;
        let lane = ctx.n_seq % LANES;
        for i in 0..n_nodes {
            let mm = m[i][v_last].to_array()[lane];
            let ee = e[i][v_last].to_array()[lane];
            let ff = f[i][v_last].to_array()[lane];
            let s = mm.max(ee).max(ff);
            if s > best_score {
                best_score = s;
                best_i = i;
                best_j = ctx.n_seq;
                best_state = if s == mm {
                    0
                } else if s == ee {
                    1
                } else {
                    2
                };
            }
        }
    }

    (m, e, f, best_score, best_i, best_j, best_state)
}

fn align_wide(
    params: &AlignmentParams,
    align_type: AlignmentType,
    sequence: &[u8],
    graph: &PoaGraph,
    ctx: &AlignCtx,
) -> Alignment {
    let (m, e, f, best_score, best_i, best_j, best_state) =
        forward_wide(params, align_type, sequence, ctx);
    let path = backtrack(
        params,
        align_type,
        &ctx.sorted_nodes,
        &ctx.node_map,
        graph,
        sequence,
        best_score,
        best_i,
        best_j,
        best_state,
        |i, j, state| {
            let v = j / LANES;
            let k = j % LANES;
            match state {
                0 => m[i][v].to_array()[k],
                1 => e[i][v].to_array()[k],
                _ => f[i][v].to_array()[k],
            }
        },
    );

    Alignment {
        score: best_score,
        path,
    }
}

#[cfg(target_arch = "x86_64")]
mod avx2 {
    use super::*;
    use std::arch::x86_64::*;

    #[inline]
    unsafe fn set1(v: i32) -> __m256i {
        _mm256_set1_epi32(v)
    }

    #[inline]
    unsafe fn load(arr: &[i32; LANES]) -> __m256i {
        _mm256_loadu_si256(arr.as_ptr() as *const __m256i)
    }

    #[inline]
    unsafe fn store_to(v: __m256i, arr: &mut [i32; LANES]) {
        _mm256_storeu_si256(arr.as_mut_ptr() as *mut __m256i, v);
    }

    /// `[prev_last, v0, ..., v6]`: shift right one i32 lane.
    #[inline]
    unsafe fn shift1(v: __m256i, prev_last: i32) -> __m256i {
        // permutevar8x32 -> [v7, v0, v1, ..., v6]; blend lane 0 with prev_last.
        let perm = _mm256_permutevar8x32_epi32(v, _mm256_setr_epi32(7, 0, 1, 2, 3, 4, 5, 6));
        _mm256_blend_epi32(perm, set1(prev_last), 0b0000_0001)
    }

    /// Prefix scan with linear penalties (log-step shift 1/2/4 lanes).
    #[inline]
    unsafe fn prefix_scan(mut a: __m256i, e: i32) -> __m256i {
        let neg_inf = set1(NEG_INF);
        for (n, excl, idx) in [
            (1, [-1, 0, 0, 0, 0, 0, 0, 0], [7, 0, 1, 2, 3, 4, 5, 6]),
            (2, [-1, -1, 0, 0, 0, 0, 0, 0], [7, 7, 0, 1, 2, 3, 4, 5]),
            (4, [-1, -1, -1, -1, 0, 0, 0, 0], [7, 7, 7, 7, 0, 1, 2, 3]),
        ] {
            // shifted[k] = a[k - n] for k >= n; excluded lanes repeat a[7]
            // (garbage, masked out below).
            let shifted = _mm256_permutevar8x32_epi32(a, load(&idx));
            let excl_mask = load(&excl);
            // Excluded lanes (k < n) become negative infinity.
            let masked = _mm256_blendv_epi8(shifted, neg_inf, excl_mask);
            a = _mm256_max_epi32(a, _mm256_add_epi32(masked, set1(n * e)));
        }
        a
    }

    #[inline]
    unsafe fn set_lane0(v: __m256i, boundary: i32) -> __m256i {
        _mm256_blend_epi32(v, set1(boundary), 0b0000_0001)
    }

    #[inline]
    unsafe fn clamp_neg_inf(v: __m256i) -> __m256i {
        let mask = _mm256_cmpgt_epi32(set1(NEG_INF_HALF), v);
        _mm256_blendv_epi8(v, set1(NEG_INF), mask)
    }

    #[target_feature(enable = "avx2")]
    #[allow(clippy::type_complexity)]
    unsafe fn forward_avx2(
        params: &AlignmentParams,
        align_type: AlignmentType,
        sequence: &[u8],
        ctx: &AlignCtx,
    ) -> (
        Vec<Vec<__m256i>>,
        Vec<Vec<__m256i>>,
        Vec<Vec<__m256i>>,
        i32,
        usize,
        usize,
        u8,
    ) {
        let m_score = params.match_score;
        let n_score = params.mismatch_score;
        let g_open = params.gap_open;
        let g_ext = params.gap_extend;
        let is_local = align_type == AlignmentType::Local;
        let is_semi = align_type == AlignmentType::SemiGlobal;

        let n_nodes = ctx.n_nodes;
        let n_vec = ctx.n_vec;
        let g_open_v = set1(g_open);
        let g_ext_v = set1(g_ext);
        let neg_inf_v = set1(NEG_INF);
        let j1_zero = load(&[
            NEG_INF, 0, NEG_INF, NEG_INF, NEG_INF, NEG_INF, NEG_INF, NEG_INF,
        ]);

        let root: Vec<__m256i> = build_root(g_open, g_ext, ctx.n_seq, n_vec)
            .iter()
            .map(|arr| load(arr))
            .collect();

        let mut m = vec![vec![neg_inf_v; n_vec]; n_nodes];
        let mut e = vec![vec![neg_inf_v; n_vec]; n_nodes];
        let mut f = vec![vec![neg_inf_v; n_vec]; n_nodes];
        let mut f0_of = vec![NEG_INF; n_nodes];

        let mut best_score = if is_local { 0 } else { NEG_INF };
        let mut best_i = 0usize;
        let mut best_j = ctx.n_seq;
        let mut best_state = 0u8;

        for i in 0..n_nodes {
            let preds = &ctx.preds_of_rank[i];
            let is_start = ctx.is_start[i];
            let profile: Vec<__m256i> = build_profile(
                sequence,
                ctx.base_of_rank[i],
                m_score,
                n_score,
                ctx.n_seq,
                n_vec,
            )
            .iter()
            .map(|arr| load(arr))
            .collect();

            let f0 = if is_local || is_semi {
                0
            } else if is_start {
                g_open
            } else {
                let mut max_prev = NEG_INF;
                for &u in preds {
                    max_prev = max_prev.max(f0_of[u]);
                }
                if max_prev > NEG_INF {
                    max_prev + g_ext
                } else {
                    NEG_INF
                }
            };
            f0_of[i] = f0;

            let mut m_last = NEG_INF;
            let mut e_last = NEG_INF;
            let mut f_last = NEG_INF;

            for v in 0..n_vec {
                let prof = profile[v];

                // M
                let mut m_vec = neg_inf_v;
                if is_start {
                    let mut root_arr = [0i32; LANES];
                    let root_last = if v == 0 {
                        NEG_INF
                    } else {
                        store_to(root[v - 1], &mut root_arr);
                        root_arr[LANES - 1]
                    };
                    m_vec = _mm256_add_epi32(shift1(root[v], root_last), prof);
                } else {
                    for &u in preds {
                        let pm = _mm256_max_epi32(_mm256_max_epi32(m[u][v], e[u][v]), f[u][v]);
                        let mut last_arr = [0i32; LANES];
                        let last = if v == 0 {
                            NEG_INF
                        } else {
                            let pv = _mm256_max_epi32(
                                _mm256_max_epi32(m[u][v - 1], e[u][v - 1]),
                                f[u][v - 1],
                            );
                            store_to(pv, &mut last_arr);
                            last_arr[LANES - 1]
                        };
                        let diag = _mm256_add_epi32(shift1(pm, last), prof);
                        m_vec = _mm256_max_epi32(m_vec, diag);
                    }
                    if (is_local || is_semi) && v == 0 {
                        m_vec = _mm256_max_epi32(m_vec, _mm256_add_epi32(j1_zero, prof));
                    }
                }
                if is_local {
                    m_vec = _mm256_max_epi32(m_vec, set1(0));
                }
                m_vec = clamp_neg_inf(m_vec);
                if v == 0 {
                    m_vec = set_lane0(m_vec, NEG_INF);
                }

                // F
                let mut f_vec = neg_inf_v;
                if !is_start {
                    for &u in preds {
                        let fm = _mm256_add_epi32(m[u][v], g_open_v);
                        let ff = _mm256_add_epi32(f[u][v], g_ext_v);
                        let fe = _mm256_add_epi32(e[u][v], g_open_v);
                        let cur = _mm256_max_epi32(_mm256_max_epi32(fm, ff), fe);
                        f_vec = _mm256_max_epi32(f_vec, cur);
                    }
                }
                if is_local {
                    let mask = _mm256_cmpgt_epi32(set1(0), f_vec);
                    f_vec = _mm256_blendv_epi8(f_vec, neg_inf_v, mask);
                }
                f_vec = clamp_neg_inf(f_vec);
                if v == 0 {
                    f_vec = set_lane0(f_vec, f0);
                }

                // E
                let base_e = _mm256_add_epi32(
                    _mm256_max_epi32(shift1(m_vec, m_last), shift1(f_vec, f_last)),
                    g_open_v,
                );
                let x = if v == 0 {
                    set_lane0(base_e, NEG_INF)
                } else {
                    let mut arr = [0i32; LANES];
                    store_to(base_e, &mut arr);
                    set_lane0(base_e, arr[0].max(e_last + g_ext))
                };
                let mut e_vec = prefix_scan(x, g_ext);
                if is_local {
                    let mask = _mm256_cmpgt_epi32(set1(0), e_vec);
                    e_vec = _mm256_blendv_epi8(e_vec, neg_inf_v, mask);
                }
                e_vec = clamp_neg_inf(e_vec);

                m[i][v] = m_vec;
                e[i][v] = e_vec;
                f[i][v] = f_vec;

                if is_local {
                    let col = _mm256_max_epi32(_mm256_max_epi32(m_vec, e_vec), f_vec);
                    let mut col_arr = [0i32; LANES];
                    store_to(col, &mut col_arr);
                    let start_k = if v == 0 { 1 } else { 0 };
                    let end_k = (ctx.n_seq + 1 - v * LANES).min(LANES);
                    if end_k > start_k {
                        let cm = col_arr[start_k..end_k].iter().copied().max().unwrap();
                        if cm >= best_score {
                            best_score = cm;
                            let mut m_arr = [0i32; LANES];
                            let mut e_arr = [0i32; LANES];
                            store_to(m_vec, &mut m_arr);
                            store_to(e_vec, &mut e_arr);
                            for k in (start_k..end_k).rev() {
                                if col_arr[k] == cm {
                                    best_i = i;
                                    best_j = v * LANES + k;
                                    best_state = if m_arr[k] == cm {
                                        0
                                    } else if e_arr[k] == cm {
                                        1
                                    } else {
                                        2
                                    };
                                    break;
                                }
                            }
                        }
                    }
                }

                let mut tail = [0i32; LANES];
                store_to(m_vec, &mut tail);
                m_last = tail[LANES - 1];
                store_to(e_vec, &mut tail);
                e_last = tail[LANES - 1];
                store_to(f_vec, &mut tail);
                f_last = tail[LANES - 1];
            }
        }

        if !is_local {
            let v_last = ctx.n_seq / LANES;
            let lane = ctx.n_seq % LANES;
            for i in 0..n_nodes {
                let mut m_arr = [0i32; LANES];
                let mut e_arr = [0i32; LANES];
                let mut f_arr = [0i32; LANES];
                store_to(m[i][v_last], &mut m_arr);
                store_to(e[i][v_last], &mut e_arr);
                store_to(f[i][v_last], &mut f_arr);
                let mm = m_arr[lane];
                let ee = e_arr[lane];
                let ff = f_arr[lane];
                let s = mm.max(ee).max(ff);
                if s > best_score {
                    best_score = s;
                    best_i = i;
                    best_j = ctx.n_seq;
                    best_state = if s == mm {
                        0
                    } else if s == ee {
                        1
                    } else {
                        2
                    };
                }
            }
        }

        (m, e, f, best_score, best_i, best_j, best_state)
    }

    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn align_avx2(
        params: &AlignmentParams,
        align_type: AlignmentType,
        sequence: &[u8],
        graph: &PoaGraph,
        ctx: &AlignCtx,
    ) -> Alignment {
        let (m, e, f, best_score, best_i, best_j, best_state) =
            forward_avx2(params, align_type, sequence, ctx);
        let path = backtrack(
            params,
            align_type,
            &ctx.sorted_nodes,
            &ctx.node_map,
            graph,
            sequence,
            best_score,
            best_i,
            best_j,
            best_state,
            |i, j, state| {
                let v = j / LANES;
                let k = j % LANES;
                let mut arr = [0i32; LANES];
                match state {
                    0 => {
                        store_to(m[i][v], &mut arr);
                        arr[k]
                    }
                    1 => {
                        store_to(e[i][v], &mut arr);
                        arr[k]
                    }
                    _ => {
                        store_to(f[i][v], &mut arr);
                        arr[k]
                    }
                }
            },
        );

        Alignment {
            score: best_score,
            path,
        }
    }
}

#[cfg(target_arch = "x86_64")]
use avx2::align_avx2;

/// Traceback shared by both SIMD paths; `get(i, j, state)` reads a matrix cell
/// (0 = M, 1 = E, 2 = F). Mirrors the scalar engine's backtracking exactly.
#[allow(clippy::too_many_arguments)]
fn backtrack<G>(
    params: &AlignmentParams,
    align_type: AlignmentType,
    sorted_nodes: &[NodeIndex],
    node_map: &HashMap<NodeIndex, usize>,
    graph: &PoaGraph,
    sequence: &[u8],
    best_score: i32,
    best_node_idx: usize,
    best_col: usize,
    best_state: u8,
    get: G,
) -> Vec<(Option<usize>, Option<NodeIndex>)>
where
    G: Fn(usize, usize, u8) -> i32,
{
    let is_local = align_type == AlignmentType::Local;
    let is_semi = align_type == AlignmentType::SemiGlobal;

    let mut path = Vec::new();
    let mut curr_i = best_node_idx;
    let mut curr_j = best_col;
    let mut curr_state = best_state;

    while curr_j > 0 || curr_i > 0 {
        let node_idx = sorted_nodes[curr_i];
        let preds: Vec<NodeIndex> = graph
            .graph
            .neighbors_directed(node_idx, petgraph::Direction::Incoming)
            .collect();
        let is_start = preds.is_empty();

        if is_local && best_score == 0 {
            break;
        }
        if is_local {
            let s = get(curr_i, curr_j, curr_state);
            if s <= 0 {
                break;
            }
        }
        if is_semi && curr_j == 0 {
            break;
        }
        if curr_j == 0 && is_start {
            break;
        }

        match curr_state {
            0 => {
                let match_score = if curr_j > 0 {
                    if sequence[curr_j - 1] == graph.graph[node_idx].base {
                        params.match_score
                    } else {
                        params.mismatch_score
                    }
                } else {
                    0
                };

                if (is_local || is_semi) && curr_j == 1 && get(curr_i, curr_j, 0) == match_score {
                    path.push((Some(curr_j - 1), Some(node_idx)));
                    curr_j -= 1;
                    break;
                }

                if is_start {
                    if curr_j > 0 {
                        path.push((Some(curr_j - 1), Some(node_idx)));
                        curr_j -= 1;
                    }
                    break;
                } else {
                    let mut found = false;
                    for &pred in &preds {
                        let u = node_map[&pred];
                        let target = get(curr_i, curr_j, 0) - match_score;

                        if get(u, curr_j - 1, 0) == target {
                            path.push((Some(curr_j - 1), Some(node_idx)));
                            curr_i = u;
                            curr_j -= 1;
                            curr_state = 0;
                            found = true;
                            break;
                        }
                        if get(u, curr_j - 1, 1) == target {
                            path.push((Some(curr_j - 1), Some(node_idx)));
                            curr_i = u;
                            curr_j -= 1;
                            curr_state = 1;
                            found = true;
                            break;
                        }
                        if get(u, curr_j - 1, 2) == target {
                            path.push((Some(curr_j - 1), Some(node_idx)));
                            curr_i = u;
                            curr_j -= 1;
                            curr_state = 2;
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        break;
                    }
                }
            }
            1 => {
                let target = get(curr_i, curr_j, 1);
                let score_e = get(curr_i, curr_j - 1, 1) + params.gap_extend;

                path.push((Some(curr_j - 1), None));

                if target == score_e {
                    curr_j -= 1;
                    curr_state = 1;
                } else {
                    let score_m = get(curr_i, curr_j - 1, 0) + params.gap_open;
                    if target == score_m {
                        curr_j -= 1;
                        curr_state = 0;
                    } else {
                        curr_j -= 1;
                        curr_state = 2;
                    }
                }
            }
            _ => {
                let mut found = false;
                for &pred in &preds {
                    let u = node_map[&pred];
                    let target = get(curr_i, curr_j, 2);
                    if get(u, curr_j, 2) + params.gap_extend == target {
                        path.push((None, Some(node_idx)));
                        curr_i = u;
                        curr_state = 2;
                        found = true;
                        break;
                    }
                    if get(u, curr_j, 0) + params.gap_open == target {
                        path.push((None, Some(node_idx)));
                        curr_i = u;
                        curr_state = 0;
                        found = true;
                        break;
                    }
                    if get(u, curr_j, 1) + params.gap_open == target {
                        path.push((None, Some(node_idx)));
                        curr_i = u;
                        curr_state = 1;
                        found = true;
                        break;
                    }
                }
                if !found {
                    if is_start {
                        path.push((None, Some(node_idx)));
                    }
                    break;
                }
            }
        }
    }

    if !is_local && !is_semi {
        while curr_j > 0 {
            path.push((Some(curr_j - 1), None));
            curr_j -= 1;
        }
    }

    path.reverse();
    path
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::poa::align::ScalarAlignmentEngine;
    use rand::{rngs::StdRng, Rng, SeedableRng};

    fn random_dag(rng: &mut StdRng, max_nodes: usize) -> PoaGraph {
        let bases = *b"ACGT";
        let mut graph = PoaGraph::new();
        let n = rng.random_range(1..=max_nodes);
        let nodes: Vec<NodeIndex> = (0..n)
            .map(|_| graph.add_node(bases[rng.random_range(0..4)]))
            .collect();
        // Edges only from lower to higher index -> guaranteed DAG.
        for i in 0..n {
            for j in i + 1..n {
                if rng.random_bool(0.35) {
                    graph.add_edge(nodes[i], nodes[j], 1);
                }
            }
        }
        graph
    }

    fn random_seq(rng: &mut StdRng, max_len: usize) -> Vec<u8> {
        let bases = *b"ACGT";
        let len = rng.random_range(0..=max_len);
        (0..len).map(|_| bases[rng.random_range(0..4)]).collect()
    }

    fn run_case(params: AlignmentParams, align_type: AlignmentType, seq: &[u8], graph: &PoaGraph) {
        let scalar = ScalarAlignmentEngine::new(params.clone(), align_type).align(seq, graph);
        let wide = align_wide(&params, align_type, seq, graph, &build_ctx(seq, graph));
        assert_eq!(
            wide.score,
            scalar.score,
            "wide score mismatch: type={align_type:?} params={params:?} seq_len={} nodes={}",
            seq.len(),
            graph.num_nodes()
        );
        assert_eq!(
            wide.path,
            scalar.path,
            "wide path mismatch: type={align_type:?} params={params:?} seq_len={} nodes={}",
            seq.len(),
            graph.num_nodes()
        );

        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx2") {
            let avx2 =
                unsafe { align_avx2(&params, align_type, seq, graph, &build_ctx(seq, graph)) };
            assert_eq!(
                avx2.score,
                scalar.score,
                "avx2 score mismatch: type={align_type:?} params={params:?} seq_len={} nodes={}",
                seq.len(),
                graph.num_nodes()
            );
            assert_eq!(
                avx2.path,
                scalar.path,
                "avx2 path mismatch: type={align_type:?} params={params:?} seq_len={} nodes={}",
                seq.len(),
                graph.num_nodes()
            );
        }
    }

    #[test]
    fn simd_matches_scalar_random() {
        let mut rng = StdRng::seed_from_u64(42);
        let param_sets = [
            AlignmentParams::default(),
            AlignmentParams {
                match_score: 2,
                mismatch_score: -3,
                gap_open: -5,
                gap_extend: -1,
            },
            AlignmentParams {
                match_score: 10,
                mismatch_score: -2,
                gap_open: -12,
                gap_extend: -8,
            },
        ];
        for &align_type in &[
            AlignmentType::Global,
            AlignmentType::Local,
            AlignmentType::SemiGlobal,
        ] {
            for params in &param_sets {
                for _ in 0..200 {
                    let graph = random_dag(&mut rng, 20);
                    let seq = random_seq(&mut rng, 40);
                    run_case(params.clone(), align_type, &seq, &graph);
                }
            }
        }
    }

    #[test]
    fn simd_matches_scalar_long_seq() {
        // Longer sequences exercise cross-vector-column diagonals and gaps.
        let mut rng = StdRng::seed_from_u64(7);
        for _ in 0..20 {
            let graph = random_dag(&mut rng, 25);
            let seq = random_seq(&mut rng, 130);
            run_case(
                AlignmentParams::default(),
                AlignmentType::Global,
                &seq,
                &graph,
            );
            run_case(
                AlignmentParams::default(),
                AlignmentType::Local,
                &seq,
                &graph,
            );
            run_case(
                AlignmentParams::default(),
                AlignmentType::SemiGlobal,
                &seq,
                &graph,
            );
        }
    }

    #[test]
    fn simd_matches_scalar_empty_seq() {
        let mut rng = StdRng::seed_from_u64(1);
        for _ in 0..10 {
            let graph = random_dag(&mut rng, 10);
            run_case(
                AlignmentParams::default(),
                AlignmentType::Global,
                &[],
                &graph,
            );
            run_case(
                AlignmentParams::default(),
                AlignmentType::Local,
                &[],
                &graph,
            );
        }
    }

    #[test]
    fn simd_matches_scalar_column_boundaries() {
        // Sequence lengths at vector-column boundaries: exact multiples of 8
        // (no padding) and 8n+1 (1 valid lane + 7 padding lanes).
        let mut rng = StdRng::seed_from_u64(99);
        let bases = *b"ACGT";
        for len in [8usize, 9, 16, 17, 64, 65, 127, 128, 129] {
            let graph = random_dag(&mut rng, 15);
            let seq: Vec<u8> = (0..len).map(|_| bases[rng.random_range(0..4)]).collect();
            for align_type in [
                AlignmentType::Global,
                AlignmentType::Local,
                AlignmentType::SemiGlobal,
            ] {
                run_case(AlignmentParams::default(), align_type, &seq, &graph);
            }
        }
    }

    #[test]
    fn simd_matches_scalar_linear_and_singleton() {
        // Linear chain (no branches) and single-node graph.
        let mut rng = StdRng::seed_from_u64(123);
        let bases = *b"ACGT";
        for n in [1usize, 2, 5, 10] {
            let mut graph = PoaGraph::new();
            let nodes: Vec<NodeIndex> = (0..n)
                .map(|_| graph.add_node(bases[rng.random_range(0..4)]))
                .collect();
            for w in nodes.windows(2) {
                graph.add_edge(w[0], w[1], 1);
            }
            let seq = random_seq(&mut rng, 30);
            for align_type in [
                AlignmentType::Global,
                AlignmentType::Local,
                AlignmentType::SemiGlobal,
            ] {
                run_case(AlignmentParams::default(), align_type, &seq, &graph);
            }
        }
    }
}
