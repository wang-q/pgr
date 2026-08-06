//! RMBlast repeat-masking recipes ported from RepeatMasker 4.2.4.

use std::io::BufRead;

/// RepeatMasker `general_search_parameters` minscore, overridable by
/// RepeatMasker's `-cutoff` (default 225). `runStage` passes it through
/// unchanged (the 7.5% discount only exists in the unused `runTestStage`).
pub const CUTOFF_DEFAULT: i32 = 225;

/// Effective `-min_raw_gapped_score`: the cutoff itself.
pub fn effective_minscore(cutoff: i32) -> i32 {
    cutoff
}

/// `minmatch` tiers selected by RepeatMasker speed options:
/// `-s` (slow) / default / `-q` (quick) / `-qq` (rush).
pub const WORD_SIZE_TIERS: [usize; 4] = [8, 9, 11, 13];

/// `-gapopen` value derived from the recipe gap_initValue (-30) and
/// ins_gap_extValue (-6): abs(-30 - -6).
const GAPOPEN: i32 = 24;

/// `-gapextend` value derived from ins_gap_extValue (-6).
const GAPEXTEND: i32 = 6;

/// `-mask_level` hard-coded by `runStage` (not the recipe value 90).
const MASK_LEVEL: i32 = 101;

/// rmblastn tab output columns (RMBlast >= 2.13, no qseq/sseq).
pub const OUTFMT: &str =
    "6 score perc_sub perc_query_gap perc_db_gap qseqid qstart qend qlen sstrand sseqid \
     sstart send slen kdiv cpg_kdiv transi transv cpg_sites";

/// RepeatMasker `runTRFStage` PERFECT parameters (young simple repeats).
pub const TRF_PERFECT_ARGS: [&str; 7] = ["2", "7", "7", "80", "10", "50", "10"];

/// RepeatMasker PERFECT-stage minimum copy number (kept when `> value`).
pub const TRF_PERFECT_MIN_COPY: f64 = 4.0;

/// RepeatMasker `runTRFStage` DIVERGED parameters (old simple repeats).
pub const TRF_DIVERGED_ARGS: [&str; 7] = ["2", "3", "5", "75", "20", "33", "7"];

/// RepeatMasker DIVERGED-stage minimum copy number (kept when `> value`).
pub const TRF_DIVERGED_MIN_COPY: f64 = 5.0;

/// The GC-keyed nucleotide matrices shipped with RepeatMasker
/// (Matrices/ncbi/nt/20p##g.matrix), embedded like the lastz presets.
pub const MATRIX_GC_NAMES: [&str; 10] = [
    "35g", "37g", "39g", "41g", "43g", "45g", "47g", "49g", "51g", "53g",
];

/// `20p35g.matrix` (RepeatMasker Matrices/ncbi/nt).
pub const MATRIX_20P35G: &str = r#"# FREQS A 0.325 C 0.175 G 0.175 T 0.325
    A   R   G   C   Y   T   K   M   S   W   N   X
A   8   1  -4 -14 -15 -17 -11  -2  -9  -4  -1 -30
R   0   1   3 -13 -15 -16  -6  -6  -5  -8  -1 -30
G  -7   1  11 -13 -14 -15  -1 -10   0 -11  -1 -30
C -15 -14 -13  11   1  -7 -10  -1   0 -11  -1 -30
Y -16 -15 -13   3   1   0  -6  -6  -5  -8  -1 -30
T -17 -15 -14  -4   1   8  -2 -11  -9  -4  -1 -30
K -12  -6  -1  -9  -6  -3  -2 -10  -5  -8  -1 -30
M  -3  -6  -9  -1  -6 -12 -10  -2  -5  -8  -1 -30
S -11  -6   0   0  -6 -11  -6  -6   0 -11  -1 -30
W  -4  -6  -9  -9  -6  -4  -6  -6  -9  -4  -1 -30
N  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
X -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30
"#;

/// `20p37g.matrix` (RepeatMasker Matrices/ncbi/nt).
pub const MATRIX_20P37G: &str = r#"# FREQS A 0.315 C 0.185 G 0.185 T 0.315
    A   R   G   C   Y   T   K   M   S   W   N   X
A   8   2  -4 -14 -16 -17 -11  -2  -9  -4  -1 -30
R   0   1   3 -14 -15 -16  -6  -7  -5  -8  -1 -30
G  -8   1  11 -14 -14 -15  -2 -11  -1 -11  -1 -30
C -15 -14 -14  11   1  -7 -10  -1  -1 -11  -1 -30
Y -16 -15 -14   3   1   0  -6  -6  -5  -8  -1 -30
T -17 -16 -14  -4   2   8  -2 -11  -9  -4  -1 -30
K -13  -7  -1  -9  -6  -3  -2 -11  -5  -8  -1 -30
M  -3  -6  -9  -1  -7 -12 -11  -2  -5  -8  -1 -30
S -11  -6  -1  -1  -6 -11  -6  -6  -1 -11  -1 -30
W  -4  -7  -9  -9  -7  -4  -7  -7  -9  -4  -1 -30
N  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
X -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30
"#;

/// `20p39g.matrix` (RepeatMasker Matrices/ncbi/nt).
pub const MATRIX_20P39G: &str = r#"# FREQS A 0.305 C 0.195 G 0.195 T 0.305
    A   R   G   C   Y   T   K   M   S   W   N   X
A   8   2  -4 -14 -16 -17 -11  -2  -9  -4  -1 -30
R   0   1   3 -14 -15 -16  -6  -7  -5  -8  -1 -30
G  -8   1  10 -14 -14 -15  -2 -11  -1 -11  -1 -30
C -15 -14 -14  10   1  -8 -11  -2  -1 -11  -1 -30
Y -16 -15 -14   3   1   0  -7  -6  -5  -8  -1 -30
T -17 -16 -14  -4   2   8  -2 -11  -9  -4  -1 -30
K -13  -7  -1  -9  -6  -3  -2 -11  -5  -8  -1 -30
M  -3  -6  -9  -1  -7 -13 -11  -2  -5  -8  -1 -30
S -11  -6  -1  -1  -6 -11  -6  -6  -1 -11  -1 -30
W  -4  -6  -9  -9  -6  -4  -6  -6  -9  -4  -1 -30
N  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
X -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30
"#;

/// `20p41g.matrix` (RepeatMasker Matrices/ncbi/nt).
pub const MATRIX_20P41G: &str = r#"# FREQS A 0.295 C 0.205 G 0.205 T 0.295
    A   R   G   C   Y   T   K   M   S   W   N   X
A   9   2  -4 -15 -16 -17 -11  -3  -9  -4  -1 -30
R   0   1   3 -14 -15 -16  -6  -7  -5  -8  -1 -30
G  -8   1  10 -14 -15 -15  -2 -11  -2 -11  -1 -30
C -15 -15 -14  10   1  -8 -11  -2  -2 -11  -1 -30
Y -16 -15 -14   3   1   0  -7  -6  -5  -8  -1 -30
T -17 -16 -15  -4   2   9  -3 -11  -9  -4  -1 -30
K -13  -7  -2  -9  -6  -3  -2 -11  -5  -8  -1 -30
M  -3  -6  -9  -2  -7 -13 -11  -2  -5  -8  -1 -30
S -11  -6  -2  -2  -6 -11  -6  -6  -2 -11  -1 -30
W  -4  -7  -9  -9  -7  -4  -7  -7  -9  -4  -1 -30
N  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
X -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30
"#;

/// `20p43g.matrix` (RepeatMasker Matrices/ncbi/nt).
pub const MATRIX_20P43G: &str = r#"# FREQS A 0.285 C 0.215 G 0.215 T 0.285
    A   R   G   C   Y   T   K   M   S   W   N   X
A   9   2  -4 -15 -16 -17 -11  -2  -9  -4  -1 -30
R   0   1   2 -15 -15 -16  -6  -7  -6  -8  -1 -30
G  -8   1  10 -15 -15 -15  -2 -11  -2 -11  -1 -30
C -15 -15 -15  10   1  -8 -11  -2  -2 -11  -1 -30
Y -16 -15 -15   2   1   0  -7  -6  -6  -7  -1 -30
T -17 -16 -15  -4   2   9  -2 -11  -9  -4  -1 -30
K -13  -7  -2  -9  -6  -2  -2 -11  -6  -8  -1 -30
M  -2  -6  -9  -2  -7 -12 -11  -2  -6  -7  -1 -30
S -11  -7  -2  -2  -7 -11  -7  -7  -2 -11  -1 -30
W  -4  -7  -9  -9  -7  -4  -7  -7  -9  -4  -1 -30
N  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
X -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30
"#;

/// `20p45g.matrix` (RepeatMasker Matrices/ncbi/nt).
pub const MATRIX_20P45G: &str = r#"# FREQS A 0.275 C 0.225 G 0.225 T 0.275
    A   R   G   C   Y   T   K   M   S   W   N   X
A   9   2  -4 -15 -16 -17 -10  -2  -9  -3  -1 -30
R   0   1   2 -15 -15 -16  -6  -7  -6  -7  -1 -30
G  -8   1  10 -15 -15 -14  -2 -11  -2 -11  -1 -30
C -14 -15 -15  10   1  -8 -11  -2  -2 -11  -1 -30
Y -16 -15 -15   2   1   0  -7  -6  -6  -7  -1 -30
T -17 -16 -15  -4   2   9  -2 -10  -9  -3  -1 -30
K -12  -7  -2  -9  -6  -2  -2 -11  -6  -7  -1 -30
M  -2  -6  -9  -2  -7 -12 -11  -2  -6  -7  -1 -30
S -11  -7  -2  -2  -7 -11  -7  -7  -2 -11  -1 -30
W  -3  -6  -9  -9  -6  -3  -6  -6  -9  -3  -1 -30
N  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
X -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30
"#;

/// `20p47g.matrix` (RepeatMasker Matrices/ncbi/nt).
pub const MATRIX_20P47G: &str = r#"# FREQS A 0.265 C 0.235 G 0.235 T 0.265
    A   R   G   C   Y   T   K   M   S   W   N   X
A   9   2  -4 -15 -15 -16 -10  -2  -9  -3  -1 -30
R   0   1   2 -15 -15 -15  -6  -7  -6  -7  -1 -30
G  -7   0   9 -15 -14 -14  -2 -11  -2 -11  -1 -30
C -14 -14 -15   9   0  -7 -11  -2  -2 -11  -1 -30
Y -15 -15 -15   2   1   0  -7  -6  -6  -7  -1 -30
T -16 -15 -15  -4   2   9  -2 -10  -9  -3  -1 -30
K -12  -7  -2  -9  -6  -2  -2 -10  -6  -7  -1 -30
M  -2  -6  -9  -2  -7 -12 -10  -2  -6  -7  -1 -30
S -11  -6  -2  -2  -6 -11  -6  -6  -2 -11  -1 -30
W  -3  -6  -9  -9  -6  -3  -6  -6  -9  -3  -1 -30
N  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
X -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30
"#;

/// `20p49g.matrix` (RepeatMasker Matrices/ncbi/nt).
pub const MATRIX_20P49G: &str = r#"# FREQS A 0.255 C 0.245 G 0.245 T 0.255
    A   R   G   C   Y   T   K   M   S   W   N   X
A   9   2  -3 -14 -15 -15  -9  -2  -9  -3  -1 -30
R   1   1   2 -15 -14 -14  -6  -6  -6  -6  -1 -30
G  -7   1   9 -15 -14 -13  -2 -11  -2 -10  -1 -30
C -13 -14 -15   9   1  -7 -11  -2  -2 -10  -1 -30
Y -14 -14 -15   2   1   1  -6  -6  -6  -6  -1 -30
T -15 -15 -14  -3   2   9  -2  -9  -9  -3  -1 -30
K -11  -7  -2  -9  -5  -2  -2 -10  -6  -6  -1 -30
M  -2  -5  -9  -2  -7 -11 -10  -2  -6  -6  -1 -30
S -10  -6  -2  -2  -6 -10  -6  -6  -2 -10  -1 -30
W  -3  -6  -9  -9  -6  -3  -6  -6  -9  -3  -1 -30
N  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
X -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30
"#;

/// `20p51g.matrix` (RepeatMasker Matrices/ncbi/nt).
pub const MATRIX_20P51G: &str = r#"# FREQS A 0.245 C 0.255 G 0.255 T 0.245
    A   R   G   C   Y   T   K   M   S   W   N   X
A   9   3  -3 -14 -14 -14  -9  -2  -9  -2  -1 -30
R   1   2   2 -14 -14 -13  -5  -6  -5  -6  -1 -30
G  -6   1   9 -14 -13 -13  -1 -10  -2  -9  -1 -30
C -13 -13 -14   9   1  -6 -10  -2  -2  -9  -1 -30
Y -14 -14 -14   2   2   1  -6  -5  -5  -6  -1 -30
T -14 -14 -14  -3   3   9  -2  -9  -9  -2  -1 -30
K -10  -6  -2  -9  -5  -1  -2 -10  -5  -6  -1 -30
M  -1  -5  -9  -2  -6 -10  -9  -2  -5  -6  -1 -30
S -10  -6  -2  -2  -6  -9  -6  -6  -2  -9  -1 -30
W  -2  -5  -9  -9  -5  -2  -5  -5  -9  -2  -1 -30
N  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
X -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30
"#;

/// `20p53g.matrix` (RepeatMasker Matrices/ncbi/nt).
pub const MATRIX_20P53G: &str = r#"# FREQS A 0.235 C 0.265 G 0.265 T 0.235
    A   R   G   C   Y   T   K   M   S   W   N   X
A   9   3  -3 -14 -14 -13  -8  -2  -8  -1  -1 -30
R   1   2   2 -14 -13 -13  -5  -6  -5  -5  -1 -30
G  -6   1   8 -14 -13 -12  -1 -10  -2  -9  -1 -30
C -12 -13 -14   8   1  -6 -10  -1  -2  -9  -1 -30
Y -13 -13 -14   2   2   1  -6  -5  -5  -5  -1 -30
T -13 -14 -14  -3   3   9  -2  -8  -8  -1  -1 -30
K -10  -6  -2  -8  -5  -1  -2  -9  -5  -5  -1 -30
M  -1  -5  -8  -2  -6 -10  -9  -2  -5  -5  -1 -30
S  -9  -6  -2  -2  -6  -9  -6  -6  -2  -9  -1 -30
W  -1  -5  -8  -8  -5  -1  -5  -5  -8  -1  -1 -30
N  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
X -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30
"#;

/// Matrix file content for a GC level (e.g. `"43g"`).
pub fn matrix_content(name: &str) -> Option<&'static str> {
    Some(match name {
        "35g" => MATRIX_20P35G,
        "37g" => MATRIX_20P37G,
        "39g" => MATRIX_20P39G,
        "41g" => MATRIX_20P41G,
        "43g" => MATRIX_20P43G,
        "45g" => MATRIX_20P45G,
        "47g" => MATRIX_20P47G,
        "49g" => MATRIX_20P49G,
        "51g" => MATRIX_20P51G,
        "53g" => MATRIX_20P53G,
        _ => return None,
    })
}
/// RepeatMasker `chooseMatrices`: map a GC percentage to a matrix level.
pub fn matrix_name_for_gc(gc_frac: i64) -> &'static str {
    if gc_frac <= 36 {
        "35g"
    } else if gc_frac <= 38 {
        "37g"
    } else if gc_frac <= 40 {
        "39g"
    } else if gc_frac <= 42 {
        "41g"
    } else if gc_frac <= 44 {
        "43g"
    } else if gc_frac <= 46 {
        "45g"
    } else if gc_frac <= 48 {
        "47g"
    } else if gc_frac <= 50 {
        "49g"
    } else if gc_frac <= 52 {
        "51g"
    } else {
        "53g"
    }
}

/// GC percentage of raw sequence bytes with RepeatMasker `getBatchAverageGC`
/// semantics: (G+C+S) / (length - XNRYMK), rounded to an integer.
pub fn gc_bytes(seq: &[u8]) -> i64 {
    let mut gc = 0i64;
    let mut ambiguous = 0i64;
    let mut len = 0i64;
    for b in seq {
        match b.to_ascii_uppercase() {
            b'G' | b'C' | b'S' => gc += 1,
            b'X' | b'N' | b'R' | b'Y' | b'M' | b'K' => ambiguous += 1,
            _ => {}
        }
        len += 1;
    }
    let gc_size = len - ambiguous;
    if gc_size == 0 {
        0
    } else {
        ((gc as f64 / gc_size as f64) * 100.0).round() as i64
    }
}

/// GC percentage of a FASTA reader (see [`gc_bytes`]).
pub fn gc_percent<R: BufRead>(reader: R) -> i64 {
    let mut seq: Vec<u8> = Vec::new();
    for line in reader.lines().map_while(Result::ok) {
        let line = line.trim();
        if line.is_empty() || line.starts_with('>') {
            continue;
        }
        seq.extend_from_slice(line.as_bytes());
    }
    gc_bytes(&seq)
}

/// Build the rmblastn command-line arguments for one query against a library
/// database. `matrix` is a GC level (`"43g"`); the file must exist in the
/// directory named by the `BLASTMAT` environment variable.
pub fn build_args(
    db: &str,
    query: &str,
    matrix: &str,
    minscore: i32,
    word_size: usize,
    out: &str,
    num_threads: usize,
) -> Vec<String> {
    vec![
        "-num_alignments".into(),
        "9999999".into(),
        "-db".into(),
        db.into(),
        "-query".into(),
        query.into(),
        "-gapopen".into(),
        GAPOPEN.to_string(),
        "-gapextend".into(),
        GAPEXTEND.to_string(),
        "-mask_level".into(),
        MASK_LEVEL.to_string(),
        "-complexity_adjust".into(),
        "-word_size".into(),
        word_size.to_string(),
        "-xdrop_ungap".into(),
        (minscore * 2).to_string(),
        "-xdrop_gap_final".into(),
        minscore.to_string(),
        "-xdrop_gap".into(),
        (minscore / 2).to_string(),
        "-min_raw_gapped_score".into(),
        minscore.to_string(),
        "-dust".into(),
        "no".into(),
        "-outfmt".into(),
        OUTFMT.into(),
        "-num_threads".into(),
        num_threads.to_string(),
        "-matrix".into(),
        format!("20p{matrix}.matrix"),
        "-out".into(),
        out.into(),
    ]
}

/// Parse one rmblastn tab row into `(qseqid, qstart, qend)` with 1-based
/// inclusive query (genome) coordinates.
pub fn parse_tab_row(line: &str) -> Option<(&str, i64, i64)> {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 7 {
        return None;
    }
    let qstart = cols[5].parse().ok()?;
    let qend = cols[6].parse().ok()?;
    Some((cols[4], qstart, qend))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn minscore_passes_through_unchanged() {
        assert_eq!(effective_minscore(225), 225);
        assert_eq!(effective_minscore(300), 300);
    }

    #[test]
    fn gc_matrix_mapping_matches_choose_matrices() {
        assert_eq!(matrix_name_for_gc(36), "35g");
        assert_eq!(matrix_name_for_gc(37), "37g");
        assert_eq!(matrix_name_for_gc(43), "43g");
        assert_eq!(matrix_name_for_gc(50), "49g");
        assert_eq!(matrix_name_for_gc(52), "51g");
        assert_eq!(matrix_name_for_gc(53), "53g");
        assert_eq!(matrix_name_for_gc(99), "53g");
    }

    #[test]
    fn gc_percent_counts_gc_over_full_length() {
        let fa = ">chr\nAACCGGTTNN\n".to_string();
        // G+C=4 over 8 non-ambiguous bases (N excluded), rounded: 50%.
        assert_eq!(gc_percent(Cursor::new(fa)), 50);
    }

    #[test]
    fn tab_row_parses_genome_coordinates() {
        let row = "123\t10.2\t0.5\t1.1\tNC_000913\t100\t250\t1000\t+\tIS1\t1\t151\t768\t12.3\t9.9\t5\t7\t6";
        assert_eq!(parse_tab_row(row), Some(("NC_000913", 100, 250)));
        assert_eq!(parse_tab_row("short"), None);
    }
}
