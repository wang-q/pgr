//! FASTA assembly statistics (N50, Nx, E-size).

/// Results of N50 and related assembly statistics.
#[derive(Debug, Clone)]
pub struct N50Stats {
    /// Number of records.
    pub record_cnt: usize,
    /// Total size (sum of all sequence lengths).
    pub total_size: usize,
    /// E-size (from GAGE).
    pub e_size: f64,
    /// Nx values, parallel to the input `opt_nx`.
    pub nx_sizes: Vec<usize>,
}

/// Calculate N50 and related statistics from sequence lengths.
///
/// `opt_nx` specifies the percentages (e.g. `[50, 90]` for N50 and N90).
/// `opt_genome` overrides the total size used for the Nx goal; when `None`,
/// the sum of `lens` is used.
pub fn calc_n50_stats(
    mut lens: Vec<usize>,
    opt_nx: &[usize],
    opt_genome: Option<usize>,
) -> N50Stats {
    lens.sort_unstable_by(|a, b| b.cmp(a));

    let record_cnt = lens.len();
    let total_size: usize = lens.iter().sum();

    // reach n_given% of total_size or genome_size
    let genome_size = opt_genome.unwrap_or(total_size);
    let goals: Vec<usize> = opt_nx
        .iter()
        .map(|el| ((*el as f64) * (genome_size as f64) / 100.0) as usize)
        .collect();

    let mut cumul_size = 0; // the cumulative size
    let mut e_size = 0.0;
    let mut nx_sizes = vec![0; goals.len()];

    for cur_size in lens {
        let prev_cumul_size = cumul_size;
        cumul_size += cur_size;

        if cumul_size > 0 {
            e_size = (prev_cumul_size as f64) / (cumul_size as f64) * e_size
                + (cur_size as f64 * cur_size as f64) / cumul_size as f64;
        }

        for (i, goal) in goals.iter().enumerate() {
            if nx_sizes[i] == 0 && cumul_size > *goal {
                nx_sizes[i] = cur_size;
            }
        }
    }

    N50Stats {
        record_cnt,
        total_size,
        e_size,
        nx_sizes,
    }
}

/// Transpose a vector of vectors.
pub fn transpose<T>(v: Vec<Vec<T>>) -> Vec<Vec<T>> {
    if v.is_empty() {
        return Vec::new();
    }
    let len = v[0].len();
    let mut iters: Vec<_> = v.into_iter().map(|n| n.into_iter()).collect();
    (0..len)
        .map(|_| {
            iters
                .iter_mut()
                .map(|n| n.next().unwrap())
                .collect::<Vec<T>>()
        })
        .collect()
}

/// Count bases in a sequence, returning `(valid_len, [A, C, G, T, N])`.
///
/// Non-standard characters (IUPAC codes, gaps) are excluded from `len` and
/// not counted in any of the five canonical bins.
pub fn count_bases(seq: &[u8]) -> (usize, [usize; 5]) {
    let mut len = 0usize;
    let mut base_cnt = [0usize; 5];

    for &el in seq {
        let nt = crate::libs::nt::to_nt(el);
        if !matches!(nt, crate::libs::nt::Nt::Invalid) {
            len += 1;
            base_cnt[nt as usize] += 1;
        }
    }

    (len, base_cnt)
}
